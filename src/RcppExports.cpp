// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/ramcmc.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// cholupdate
arma::mat cholupdate(arma::mat L, arma::vec u);
static SEXP ramcmc_cholupdate_try(SEXP LSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(cholupdate(L, u));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP ramcmc_cholupdate(SEXP LSEXP, SEXP uSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(ramcmc_cholupdate_try(LSEXP, uSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// choldowndate
arma::mat choldowndate(arma::mat L, arma::vec u);
static SEXP ramcmc_choldowndate_try(SEXP LSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(choldowndate(L, u));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP ramcmc_choldowndate(SEXP LSEXP, SEXP uSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(ramcmc_choldowndate_try(LSEXP, uSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// adjust_S
arma::mat adjust_S(arma::mat S, arma::vec u, double current, double target, unsigned int n, double gamma);
static SEXP ramcmc_adjust_S_try(SEXP SSEXP, SEXP uSEXP, SEXP currentSEXP, SEXP targetSEXP, SEXP nSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type current(currentSEXP);
    Rcpp::traits::input_parameter< double >::type target(targetSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(adjust_S(S, u, current, target, n, gamma));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP ramcmc_adjust_S(SEXP SSEXP, SEXP uSEXP, SEXP currentSEXP, SEXP targetSEXP, SEXP nSEXP, SEXP gammaSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(ramcmc_adjust_S_try(SSEXP, uSEXP, currentSEXP, targetSEXP, nSEXP, gammaSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int ramcmc_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("arma::mat(*cholupdate)(arma::mat,arma::vec)");
        signatures.insert("arma::mat(*choldowndate)(arma::mat,arma::vec)");
        signatures.insert("arma::mat(*adjust_S)(arma::mat,arma::vec,double,double,unsigned int,double)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP ramcmc_RcppExport_registerCCallable() { 
    R_RegisterCCallable("ramcmc", "ramcmc_cholupdate", (DL_FUNC)ramcmc_cholupdate_try);
    R_RegisterCCallable("ramcmc", "ramcmc_choldowndate", (DL_FUNC)ramcmc_choldowndate_try);
    R_RegisterCCallable("ramcmc", "ramcmc_adjust_S", (DL_FUNC)ramcmc_adjust_S_try);
    R_RegisterCCallable("ramcmc", "ramcmc_RcppExport_validate", (DL_FUNC)ramcmc_RcppExport_validate);
    return R_NilValue;
}
