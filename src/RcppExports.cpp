// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/ramcmc.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// chol_updateR
arma::mat chol_updateR(arma::mat L, arma::vec u);
RcppExport SEXP ramcmc_chol_updateR(SEXP LSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(chol_updateR(L, u));
    return rcpp_result_gen;
END_RCPP
}
// chol_downdateR
arma::mat chol_downdateR(arma::mat L, arma::vec u);
RcppExport SEXP ramcmc_chol_downdateR(SEXP LSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(chol_downdateR(L, u));
    return rcpp_result_gen;
END_RCPP
}
// adapt_SR
arma::mat adapt_SR(arma::mat S, arma::vec u, double current, double target, unsigned int n, double gamma);
RcppExport SEXP ramcmc_adapt_SR(SEXP SSEXP, SEXP uSEXP, SEXP currentSEXP, SEXP targetSEXP, SEXP nSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type current(currentSEXP);
    Rcpp::traits::input_parameter< double >::type target(targetSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(adapt_SR(S, u, current, target, n, gamma));
    return rcpp_result_gen;
END_RCPP
}
