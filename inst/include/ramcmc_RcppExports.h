// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_ramcmc_RCPPEXPORTS_H_GEN_
#define RCPP_ramcmc_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace ramcmc {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("ramcmc", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("ramcmc", "ramcmc_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in ramcmc");
            }
        }
    }

    inline arma::mat cholupdate(arma::mat L, arma::vec u) {
        typedef SEXP(*Ptr_cholupdate)(SEXP,SEXP);
        static Ptr_cholupdate p_cholupdate = NULL;
        if (p_cholupdate == NULL) {
            validateSignature("arma::mat(*cholupdate)(arma::mat,arma::vec)");
            p_cholupdate = (Ptr_cholupdate)R_GetCCallable("ramcmc", "ramcmc_cholupdate");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cholupdate(Rcpp::wrap(L), Rcpp::wrap(u));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat choldowndate(arma::mat L, arma::vec u) {
        typedef SEXP(*Ptr_choldowndate)(SEXP,SEXP);
        static Ptr_choldowndate p_choldowndate = NULL;
        if (p_choldowndate == NULL) {
            validateSignature("arma::mat(*choldowndate)(arma::mat,arma::vec)");
            p_choldowndate = (Ptr_choldowndate)R_GetCCallable("ramcmc", "ramcmc_choldowndate");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_choldowndate(Rcpp::wrap(L), Rcpp::wrap(u));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat adjust_S(arma::mat S, arma::vec u, double current, double target, unsigned int n, double gamma) {
        typedef SEXP(*Ptr_adjust_S)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_adjust_S p_adjust_S = NULL;
        if (p_adjust_S == NULL) {
            validateSignature("arma::mat(*adjust_S)(arma::mat,arma::vec,double,double,unsigned int,double)");
            p_adjust_S = (Ptr_adjust_S)R_GetCCallable("ramcmc", "ramcmc_adjust_S");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_adjust_S(Rcpp::wrap(S), Rcpp::wrap(u), Rcpp::wrap(current), Rcpp::wrap(target), Rcpp::wrap(n), Rcpp::wrap(gamma));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

}

#endif // RCPP_ramcmc_RCPPEXPORTS_H_GEN_
