#include <RcppArmadillo.h>
#include "../inst/include/ramcmc.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Cholesky update
// Given the lower triangular matrix L obtained from the Cholesky decomposition of A,
// updates L such that it corresponds to the decomposition of A + u*u'.
//
// [[Rcpp::export]]
arma::mat chol_updateR(arma::mat L, arma::vec u) {
  ramcmc::chol_update(L, u);
  return L;
}


// Cholesky downdate
// Given the lower triangular matrix L obtained from the Cholesky decomposition of A,
// updates L such that it corresponds to the decomposition of A - u*u'.
//
// NOTE: The function does not check that the downdating produces a positive definite matrix!
//
// [[Rcpp::export]]
arma::mat chol_downdateR(arma::mat L, arma::vec u) {
  ramcmc::chol_downdate(L, u);
  return L;
}

// Update the Cholesky factor of the covariance matrix of the proposal distribution
// [[Rcpp::export]]
arma::mat adapt_SR(arma::mat S, arma::vec u, double current, double target, unsigned int n, double gamma) {
 ramcmc::adapt_S(S, u, current, target, n, gamma);
 return S;
}
