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
std::vector<int> test(int ndata, int nsample){
  std::vector<int> i_sampled;
  std::vector<int> i_unsampled;
  std::vector<int> idx(ndata);
  std::iota(idx.begin(),idx.end(),0);
  
  i_sampled = as<std::vector<int>>(Rcpp::sample(ndata,nsample,false));
  i_unsampled = vector_diff(idx,i_sampled);
  return(i_unsampled);
}

