//Kiri Daust, 2020

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

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
  NumericMatrix rmat(nc, nc);
  
  vector<asset_info> info(nc);
  for (int c = 0; c < nc; c++)
    info[c] = compute_asset_info(mat, c, rstart, rend);
  
  for (int c1 = 0; c1 < nc; c1++) {
    for (int c2 = 0; c2 < c1; c2++) {
      double sXY = 0;
      
      for (int r = rstart; r < rend; r++)
        sXY += mat(r, c1) * mat(r, c2);
      
      rmat(c1, c2) = (nperiod * sXY - info[c1].sum * info[c2].sum) / (info[c1].stdev * info[c2].stdev);
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
  int objRes;
  NumericVector obj_cont_res;
};

//objective function
// [[Rcpp::export]]
objResult obj_fn(mat x, NumericMatrix strata, NumericMatrix cor_full, int eta = 1){
  int num_vars = x.ncol();
  int num_obs = strata.nrow();
  NumericVector hist_cnt;
  IntegerMatrix hist_out(num_obs, num_vars);
  NumericVector data;
  NumericVector strata_curr;
  NumericVector hist_temp;
  NumericVector obj_cont;
  NumericMatrix t2;
  
  for(int i = 0; i < num_vars; i++){
    data = x(_,i);
    strata_curr = strata(_,i);
    hist_out(_,i) = hist(data,strata_curr);
  }
  
  NumericVector temp = wrap(abs(hist_out - eta));
  temp.attr("dim") = Dimension(num_obs,num_vars);
  t2 = as<NumericMatrix>(temp);
  obj_cont = rowSums(t2);
  //Rcout << "Vectout : " << obj_cont << "\n";
  NumericMatrix cor_new = c_cor(x);
  double obj_cor = sum(abs(cor_full - cor_new));
  double objFinal = sum(obj_cont) + obj_cor*2;
  struct objRes out = {objFinal, obj_cont};
  return(out);
}


//Main Function
// [[Rcpp::export]]
//[[Rcpp::depends(RcppArmadillo)]]
List CppLHS(NumericMatrix x, int ndata, int nsample, bool cost_mode, int iter, 
                        double temperature, NumericVector cost, NumericMatrix strata, NumericMatrix cor_full){
  double prev_obj;
  double obj;
  double delta_obj;
  double metropolis;
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
  int i_worse;
  IntegerVector idx2;
  NumericVector temp;
  
  i_sampled = as<std::vector<int>>(Rcpp::sample(ndata,nsample,false));
  i_unsampled = vector_diff(idx,i_sampled);
  x_curr = xA.rows(as<arma::uvec>(i_sampled)); // is this efficient?
  struct objRes res = obj_fun(x_curr,strata,cor_full); //test this
  obj = res.objRes;
  delta_cont = res.obj_cont_res;
  
  if(cost_mode){
    idx2 = wrap(i_sampled);
    temp = cost[idx2];
    opCost = sum(temp);
  }
  
  NumericVector obj_values(iter);
  
  for(int i = 0; i < iter; i++){
    prev_obj = obj;
    i_sampled_prev = i_sampled;
    i_unsampled_prev = i_unsampled;
    delta_cont_prev = delta_cont;
    
    if(cost_mode){
      prev_opCost = opCost;
    }
    
    if(Rcpp::runif(1,0,1)[0] < 0.5){
      idx_removed = as<std::vector<int>>(Rcpp::sample(i_sampled.size(), 1, false));
      spl_removed = i_sampled[idx_removed[0]];
      i_sampled.erase(i_sampled.begin() + idx_removed[0]);
      idx_added = as<std::vector<int>>(Rcpp::sample(i_unsampled.size(), 1, false));
      i_sampled.push_back(i_unsampled[idx_added[0]]);
      i_unsampled.erase(i_unsampled.begin()+idx_added[0]);
      i_unsampled.push_back(spl_removed);
      //now ready for data
    }else{
      i_worse = *std::max_element(delta_cont.begin(),delta_cont.end());
      spl_removed = i_sampled[i_worse];
      idx_added = as<std::vector<int>>(Rcpp::sample(i_unsampled.size(), 1, false));
      i_sampled.erase(i_sampled.begin()+i_worse);
      i_sampled.push_back(i_unsampled[idx_added[0]]);
      i_unsampled.erase(i_unsampled.begin()+idx_added[0]);
      i_unsampled.push_back(spl_removed);
      //ready for data
    }
    x_curr = xA.rows(as<arma::uvec>(i_sampled)); // is this efficient?
    struct objRes res = obj_fun(x_curr,strata,cor_full); //test this
    obj = res.objRes;
    delta_cont = res.obj_cont_res;
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
      i_sampled = i_sampled_prev;
      i_unsampled = i_unsampled_prev;
      obj = prev_obj;
      delta_cont = delta_cont_prev;
      if(cost_mode){
        opCost = prev_opCost;
      }
    }
    obj_values[i] = obj;
    if(i % length.cycle == 0){
      temperature = temperature*tdecrease;
    }
  }
  x_curr = xA.rows(as<arma::uvec>(i_sampled));
  return List::create(_["sampled_data"] = x_curr,
                      _["obj"] = obj_values,
                      _["final_obj"] = delta_cont);
}


