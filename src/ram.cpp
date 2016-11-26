#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

// Cholesky update
// Given the lower triangular matrix L obtained from the Cholesky decomposition of A,
// updates L such that it corresponds to the decomposition of A + u*u'.
// [[Rcpp::export]]
arma::mat cholupdate(arma::mat L, arma::vec u) {
  unsigned int n = u.n_elem - 1;
  for (arma::uword i = 0; i < n; i++) {
    double r = sqrt(L(i,i) * L(i,i) + u(i) * u(i));
    double c = r / L(i, i);
    double s = u(i) / L(i, i);
    L(i, i) = r;
    L(arma::span(i + 1, n), i) =
      (L(arma::span(i + 1, n), i) + s * u.rows(i + 1, n)) / c;
    u.rows(i + 1, n) = c * u.rows(i + 1, n) -
      s * L(arma::span(i + 1, n), i);
  }
  L(n, n) = sqrt(L(n, n) * L(n, n) + u(n) * u(n));
  return L;
}


// Cholesky downdate
// Given the lower triangular matrix L obtained from the Cholesky decomposition of A,
// updates L such that it corresponds to the decomposition of A - u*u'.
//
// NOTE: The function does not check that the downdating produces a positive definite matrix!
//       see checks on adjust_S.
//
// [[Rcpp::export]]
arma::mat choldowndate(arma::mat L, arma::vec u) {
  unsigned int n = u.n_elem - 1;
  for (arma::uword i = 0; i < n; i++) {
    double r = sqrt(L(i,i) * L(i,i) - u(i) * u(i));
    double c = r / L(i, i);
    double s = u(i) / L(i, i);
    L(i, i) = r;
    L(arma::span(i + 1, n), i) =
      (L(arma::span(i + 1, n), i) - s * u.rows(i + 1, n)) / c;
    u.rows(i + 1, n) = c * u.rows(i + 1, n) -
      s * L(arma::span(i + 1, n), i);
  }
  L(n, n) = sqrt(L(n, n) * L(n, n) - u(n) * u(n));
  return L;
}

// pass by reference
// I don't want to export this to R side but apparently I must
// in order to it to be included in inst/include...
// [[Rcpp::export]]
void adjust_S(arma::mat& S, arma::vec& u, double current, double target, unsigned int n, double gamma) {

  double change = current - target;
  u = S * u / arma::norm(u) * sqrt(std::min(1.0, u.n_elem * pow(n, -gamma)) *
    std::abs(change));

  if(change > 0.0) {
    S = cholupdate(S, u);
  } else {
    //downdate S unless numerical problems occur
    //do nothing in case of numerical issues
    arma::mat Stmp = choldowndate(S, u);
    //should stop here in case of problems
    if(Stmp.is_finite()){
      //check diagonal
      arma::uvec cond = arma::find(arma::diagvec(Stmp) < 0);
      if (cond.n_elem == 0) {
        S = Stmp;
      }
    }
  }
}
