#include <Rcpp.h>
using namespace Rcpp;

std::vector<int> vector_diff( const std::vector<int>& model, const std::vector<int>& pattern ){
  std::set<int> s_model( model.begin(), model.end() );
  std::set<int> s_pattern( pattern.begin(), pattern.end() );
  std::vector<int> result;
  
  std::set_difference( s_model.begin(), s_model.end(), s_pattern.begin(), s_pattern.end(),
                       std::back_inserter( result ) );
  
  return result;
}

// [[Rcpp::export]]
NumericVector clhs_iter(IntegerVector i_data, int nsample, bool cost_mode, int iter, int ndata){
  double prev_obj;
  double prev_opCost;
  double opCost;
  IntegerVector prev_sampled;
  IntegerVector prev_unsampled;
  NumericVector prev_contObj;
  
  std::vector<int> idx_removed;
  int spl_removed;
  std::vector<int> idx_added;
  std::vector<int> i_sampled;
  std::vector<int> i_unsampled;
  std::vector<int> idx(ndata);
  std::iota(idx.begin(),idx.end(),0);
  
  i_sampled = as<std::vector<int>>(Rcpp::sample(ndata,nsample,false));
  i_unsampled = vector_diff(idx,i_sampled);
  //need to sample previous
  
  for(int i = 0; i < iter; i++){
    //store previous
    if(cost_mode){
      prev_opCost = opCost;
    }
    
    if(Rcpp::runif(1,0,1)[0] < 0.5){
      idx_removed = as<std::vector<int>>(Rcpp::sample(i_sampled.size(), 1, false));
      spl_removed = i_sampled[idx_removed[0]];
      i_sampled.erase(idx_removed[0]);
    }
  }
}


