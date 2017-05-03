// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List quick_lm(const arma::vec & y, const arma::mat & X) {

  arma::colvec coef = arma::solve(X, y);

  return List::create(Named("coefficients") = coef);
}
