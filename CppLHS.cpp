/*
 * Kiri Daust, August 2020
 * C++ Version of clhs by Pierre Roudier
 * 
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

//difference in vectors
std::vector<int> vector_diff( const std::vector<int>& model, const std::vector<int>& pattern ){
  std::set<int> s_model( model.begin(), model.end() );
  std::set<int> s_pattern( pattern.begin(), pattern.end() );
  std::vector<int> result;
  
  std::set_difference( s_model.begin(), s_model.end(), s_pattern.begin(), s_pattern.end(),
                       std::back_inserter( result ) );
  
  return result;
}

struct asset_info {
  double sum, sum2, stdev;
};

//[correlation matrix](http://en.wikipedia.org/wiki/Correlation_and_dependence).
// n,sX,sY,sXY,sX2,sY2
// cor = ( n * sXY - sX * sY ) / ( sqrt(n * sX2 - sX^2) * sqrt(n * sY2 - sY^2) )
inline asset_info compute_asset_info(const NumericMatrix& mat, 
                                     const int icol, const int rstart, const int rend) {
  double sum, sum2;
  sum = sum2 = 0;
  
  for (int r = rstart; r < rend; r++) {
    double d = mat(r, icol);
    sum += d;
    sum2 += pow(d,2);
  }
  
  asset_info res;
  res.sum = sum;
  res.sum2 = sum2;
  res.stdev = sqrt((rend-rstart) * sum2 - pow(sum, 2));
  return res;
}

inline NumericMatrix c_cor_helper(const NumericMatrix& mat, const int rstart, const int rend) {
  int nc = mat.ncol();
  int nperiod = rend - rstart;
  double tempCor;
  NumericMatrix rmat(nc, nc);
  
  vector<asset_info> info(nc);
  for (int c = 0; c < nc; c++)
    info[c] = compute_asset_info(mat, c, rstart, rend);
  
  for (int c1 = 0; c1 < nc; c1++) {
    for (int c2 = 0; c2 < c1; c2++) {
      double sXY = 0;
      
      for (int r = rstart; r < rend; r++)
        sXY += mat(r, c1) * mat(r, c2);
      
      tempCor = (nperiod * sXY - info[c1].sum * info[c2].sum) / (info[c1].stdev * info[c2].stdev);
      if(isnan(tempCor)){
        tempCor = 0;
      }
      rmat(c1, c2) = tempCor;
    }
  }
  
  return rmat;
}

// [[Rcpp::export]]
NumericMatrix c_cor(NumericMatrix mat) {
  return c_cor_helper(mat, 0, mat.nrow());
}

//bincount
IntegerVector hist(NumericVector x, NumericVector breaks){ //based on C_bincount from graphics package
    int n = x.length();
    int nb = breaks.length();
    int nb1 = nb-1;
    int i,lo,hi,newVal;
    bool right = true,include_border = true;
    //R_xlen_t i, lo, hi, nb1 = nb - 1, new;
    
    IntegerVector counts(nb1);
    
    for(i = 0 ; i < n ; i++){
      lo = 0;
      hi = nb1;
      if(breaks[lo] <= x[i] &&
         (x[i] < breaks[hi] || (x[i] == breaks[hi]))) {
        while(hi-lo >= 2) {
          newVal = (hi+lo)/2;
          if(x[i] > breaks[newVal]){
            lo = newVal;
          }else{
            hi = newVal;
          }
        }
        counts[lo]++;
      }
    }
    return(counts);
  }

struct objResult {
  double objRes;
  std::vector<double> obj_cont_res;
};

//objective function
objResult obj_fn(arma::mat x, NumericMatrix strata, NumericMatrix cor_full, int eta = 1){
  int num_vars = x.n_cols;
  int num_obs = strata.nrow();
  NumericVector hist_cnt;
  IntegerMatrix hist_out(num_obs - 1, num_vars);
  NumericVector data;
  NumericVector strata_curr;
  NumericVector hist_temp;
  NumericVector obj_cont;
  NumericMatrix t2;
  std::vector<double> obj_cont2;
  
  for(int i = 0; i < num_vars; i++){
    data = wrap(x.col(i));
    strata_curr = strata(_,i);
    hist_out(_,i) = hist(data,strata_curr);
  }
  
  //Rcout << "Hist values \n" << hist_out << "\n";
  NumericVector temp = wrap(abs(hist_out - eta));
  temp.attr("dim") = Dimension(num_obs-1,num_vars);
  t2 = as<NumericMatrix>(temp);
  //Rcout << "New values \n" << t2 << "\n";
  obj_cont = rowSums(t2);
  obj_cont2 = as<std::vector<double>>(obj_cont);
  //Rcout << "RowSums : " << obj_cont << "\n";
  NumericMatrix cor_new = c_cor(wrap(x));
  //Rcout << "cormat " << cor_new << "\n";
  double obj_cor = sum(abs(cor_full - cor_new));
  double objFinal = sum(obj_cont) + obj_cor*2;
  struct objResult out = {objFinal, obj_cont2};
  return(out);
}


//Main Function

// [[Rcpp::export]]
List CppLHS(NumericMatrix x, NumericVector cost, NumericMatrix strata, int nsample, bool cost_mode, int iter, 
                        double temperature = 1, double tdecrease = 0.95, int length_cycle = 8){
  
  int ndata = x.nrow();
  double prev_obj;
  double obj;
  double delta_obj;
  double metropolis = 1.0;
  double prev_opCost;
  double opCost;
  double metropolis_cost;
  double delta_cost;
  IntegerVector prev_sampled;
  IntegerVector prev_unsampled;
  NumericVector prev_contObj;
  arma::mat x_curr;
  arma::mat xA = as<arma::mat>(x);
  NumericMatrix cor_mat = c_cor(x);
  
  std::vector<double> delta_cont;
  std::vector<double> delta_cont_prev;
  std::vector<int> idx_removed;
  int spl_removed;
  std::vector<int> idx_added;
  std::vector<int> i_sampled;
  std::vector<int> i_sampled_prev;
  std::vector<int> i_unsampled;
  std::vector<int> i_unsampled_prev;
  std::vector<int> idx(ndata);
  std::iota(idx.begin(),idx.end(),0);
  std::vector<double>::iterator it_worse;
  int i_worse;
  
  IntegerVector idx2;
  NumericVector temp;
  
  i_sampled = as<std::vector<int>>(Rcpp::sample(ndata-1,nsample,false));
  i_unsampled = vector_diff(idx,i_sampled);
  arma::uvec arm_isamp = arma::conv_to<arma::uvec>::from(i_sampled);
  x_curr = xA.rows(arm_isamp); // is this efficient?
  //Rcout << "Finish initial sample \n";
  NumericMatrix cor_full = c_cor(x);
  //Rcout << "cormat " << cor_full << "\n";
  struct objResult res = obj_fn(x_curr,strata,cor_full);
  obj = res.objRes;
  delta_cont = res.obj_cont_res;
  //Rcout << "Finished function; obj = "<< obj << "\n";
  
  if(cost_mode){
    idx2 = wrap(i_sampled);
    temp = cost[idx2];
    opCost = sum(temp);
  }
  
  NumericVector obj_values(iter);
  
  for(int i = 0; i < iter; i++){
    //Rcout << i << " ";
    if(i % 50 == 0)
      Rcpp::checkUserInterrupt();
    prev_obj = obj;
    i_sampled_prev = i_sampled;
    i_unsampled_prev = i_unsampled;
    delta_cont_prev = delta_cont;
    
    if(cost_mode){
      prev_opCost = opCost;
    }
    
    if(Rcpp::runif(1,0,1)[0] < 0.4){
      //Rcout << "In Simple Swap \n";
      idx_removed = as<std::vector<int>>(Rcpp::sample(i_sampled.size()-1, 1, false));
      spl_removed = i_sampled[idx_removed[0]];
      i_sampled.erase(i_sampled.begin() + idx_removed[0]);
      idx_added = as<std::vector<int>>(Rcpp::sample(i_unsampled.size()-1, 1, false));
      i_sampled.push_back(i_unsampled[idx_added[0]]);
      i_unsampled.erase(i_unsampled.begin()+idx_added[0]);
      i_unsampled.push_back(spl_removed);
      //now ready for data
    }else{
      //Rcout << "In remove worst \n";
      it_worse = std::max_element(delta_cont.begin(),delta_cont.end());
      i_worse = std::distance(delta_cont.begin(), it_worse);
      spl_removed = i_sampled[i_worse];
      //Rcout << "spl_removed " << spl_removed << "\n";
      idx_added = as<std::vector<int>>(Rcpp::sample(i_unsampled.size()-1, 1, false));
      i_sampled.erase(i_sampled.begin()+i_worse);
      i_sampled.push_back(i_unsampled[idx_added[0]]);
      i_unsampled.erase(i_unsampled.begin()+idx_added[0]);
      i_unsampled.push_back(spl_removed);
      //ready for data
    }
    //Rcout << "sending to function \n";
    arm_isamp = arma::conv_to<arma::uvec>::from(i_sampled);
    x_curr = xA.rows(arm_isamp); // is this efficient?
    res = obj_fn(x_curr,strata,cor_full); //test this
    obj = res.objRes;
    delta_cont = res.obj_cont_res;
    NumericVector dcTemp = wrap(delta_cont);
    //Rcout << "New delta obj cont is " << dcTemp << "\n";
    //update variables
    delta_obj = obj - prev_obj;
    metropolis = exp(-1*delta_obj/temperature);
    
    if(cost_mode){
      idx2 = wrap(i_sampled);
      temp = cost[idx2];
      opCost = sum(temp);
      delta_cost = opCost - prev_opCost;
      metropolis_cost = exp(-1*delta_cost/temperature);
    }else{
      metropolis_cost = R_PosInf;
    }
    
    //Revert Change
    if(delta_obj > 0 && runif(1,0,1)[0] >= metropolis || runif(1,0,1)[1] >= metropolis_cost){
      //Rcout << "In revert change \n";
      i_sampled = i_sampled_prev;
      i_unsampled = i_unsampled_prev;
      obj = prev_obj;
      delta_cont = delta_cont_prev;
      if(cost_mode){
        opCost = prev_opCost;
      }
    }
    obj_values[i] = obj;
    if(i % length_cycle == 0){
      temperature = temperature*tdecrease;
    }
  }
  arm_isamp = arma::conv_to<arma::uvec>::from(i_sampled);
  x_curr = xA.rows(arm_isamp);
  Rcout << "vroom vroom \n";
  return List::create(_["sampled_data"] = x_curr,
                      _["indeces"] = i_sampled,
                      _["obj"] = obj_values,
                      _["final_obj"] = delta_cont);
}