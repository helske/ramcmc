#include <RcppArmadillo.h>
#include "../inst/include/ramcmc.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Cholesky update
// Given the lower triangular matrix L obtained from the Cholesky decomposition of A,
// updates L such that it corresponds to the decomposition of A + u*u'.
//
// [[Rcpp::export]]
arma::mat chol_updateR(arma::mat& L, arma::vec& u) {
  return ramcmc::chol_update(L, u);
}


// Cholesky downdate
// Given the lower triangular matrix L obtained from the Cholesky decomposition of A,
// updates L such that it corresponds to the decomposition of A - u*u'.
//
// NOTE: The function does not check that the downdating produces a positive definite matrix!
//       see checks on adjust_L.
//
// [[Rcpp::export]]
arma::mat chol_downdateR(arma::mat& L, arma::vec& u) {
  return ramcmc::chol_downdate(L, u);
}


// Update the Cholesky factor of the covariance matrix of the proposal distribution
// Note that pass-by-reference, so not good for calling from R (see wrapper below)
// [[Rcpp::export]]
arma::mat adapt_LR(arma::mat L, arma::vec u, double current, double target, unsigned int n, double gamma) {
 ramcmc::adapt_L(L, u, current, target, n, gamma);
 return L;
}
