#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

std::vector<int> vector_diff( const std::vector<int>& model, const std::vector<int>& pattern ){
  std::set<int> s_model( model.begin(), model.end() );
  std::set<int> s_pattern( pattern.begin(), pattern.end() );
  std::vector<int> result;
  
  std::set_difference( s_model.begin(), s_model.end(), s_pattern.begin(), s_pattern.end(),
                       std::back_inserter( result ) );
  
  return result;
}

// [[Rcpp::export]]
CharacterVector testList(NumericMatrix x){
  List L(5);
  for(int i = 0; i < 5; i++){
    L[i] = Rcpp::table(x(_,i));
  }
  NumericVector y = L[1];
  CharacterVector z = y.attr("names");
  return(z);
}

// [[Rcpp::export]]
arma::mat test2(arma::mat x, arma::mat include){
  x = join_vert(x,include);
  return(x);
}

// // [[Rcpp::export]]
// int test(int ndata, int nsample){
//   std::vector<int> idx_removed;
//   int spl_removed;
//   std::vector<int> idx_added;
//   std::vector<int> i_sampled;
//   std::vector<int> i_unsampled;
//   std::vector<int> idx(ndata);
//   std::iota(idx.begin(),idx.end(),0);
//   
//   i_sampled = as<std::vector<int>>(Rcpp::sample(ndata,nsample,false));
//   i_unsampled = vector_diff(idx,i_sampled);
//   
//   idx_removed = as<std::vector<int>>(Rcpp::sample(i_sampled.size(), 1, false));
//   spl_removed = i_sampled[idx_removed[0]];
//   i_sampled.erase(i_sampled.begin() + idx_removed[0]);
//   idx_added = as<std::vector<int>>(Rcpp::sample(i_unsampled.size(), 1, false));
//   i_sampled.push_back(i_unsampled[idx_added[0]]);
//   int maxInd = *std::max_element(i_sampled.begin(),i_sampled.end());
//   return(maxInd);
// }


// [[Rcpp::export]]
arma::mat test1(NumericMatrix mat1, IntegerVector ind){
  arma::mat mat2 = as<arma::mat>(mat1);
  arma::uvec ind2 = as<arma::uvec>(ind);
  arma::mat mat3 = mat2.rows(ind2);
  return(mat3);
}


// double test2(int ndata, int nsample, NumericVector cost){
//   std::vector<int> idx_removed;
//   int spl_removed;
//   std::vector<int> idx_added;
//   std::vector<int> i_sampled;
//   std::vector<int> i_unsampled;
//   std::vector<int> idx(ndata);
//   std::iota(idx.begin(),idx.end(),0);
//   IntegerVector idx2;
//   NumericVector temp;
//   double out;
//   
//   i_sampled = as<std::vector<int>>(Rcpp::sample(ndata,nsample,false));
//   i_unsampled = vector_diff(idx,i_sampled);
//   idx2 = wrap(i_sampled);
//   temp = cost[idx2];
//   out = sum(temp);
//   return(out);
// }

// [[Rcpp::export]]
IntegerVector C_bincount(NumericVector x, NumericVector breaks){
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

