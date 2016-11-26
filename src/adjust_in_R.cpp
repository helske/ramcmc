#include <RcppArmadillo.h>
#include "chol.h"

// pass by value
// use this on R side so there is no side-effects
// [[Rcpp::export]]
arma::mat adjust_S_copy(arma::mat S, arma::vec u, double current, double target, unsigned int n, double gamma) {

  double change = current - target;
  u = S * u / arma::norm(u) * sqrt(std::min(1.0, u.n_elem * pow(n, -gamma)) *
    std::abs(change));

  if(change > 0.0) {
    S = cholupdate(S, u);

  } else {
    //downdate S unless numerical problems occur
    //do nothing in case of numerical issues
    arma::mat Stmp = choldowndate(S, u);
    //should stop here
    if(Stmp.is_finite()){
      //make additional check
      arma::uvec cond = arma::find(arma::diagvec(Stmp) < 0);
      if (cond.n_elem == 0) {
        S = Stmp;
      }
    }
  }

  return S;
}
