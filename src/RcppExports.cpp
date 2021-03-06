// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sigma_estim_lwident_cpp
List sigma_estim_lwident_cpp(arma::mat data, double shrink_int, bool zeromean_log);
RcppExport SEXP _CovEstim_sigma_estim_lwident_cpp(SEXP dataSEXP, SEXP shrink_intSEXP, SEXP zeromean_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type shrink_int(shrink_intSEXP);
    Rcpp::traits::input_parameter< bool >::type zeromean_log(zeromean_logSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_estim_lwident_cpp(data, shrink_int, zeromean_log));
    return rcpp_result_gen;
END_RCPP
}
// sigma_estim_lwone_cpp
List sigma_estim_lwone_cpp(arma::mat data, double shrink_int, bool zeromean_log);
RcppExport SEXP _CovEstim_sigma_estim_lwone_cpp(SEXP dataSEXP, SEXP shrink_intSEXP, SEXP zeromean_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type shrink_int(shrink_intSEXP);
    Rcpp::traits::input_parameter< bool >::type zeromean_log(zeromean_logSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_estim_lwone_cpp(data, shrink_int, zeromean_log));
    return rcpp_result_gen;
END_RCPP
}
// sigma_estim_lwcc_cpp
List sigma_estim_lwcc_cpp(arma::mat data, double shrink_int, bool zeromean_log);
RcppExport SEXP _CovEstim_sigma_estim_lwcc_cpp(SEXP dataSEXP, SEXP shrink_intSEXP, SEXP zeromean_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type shrink_int(shrink_intSEXP);
    Rcpp::traits::input_parameter< bool >::type zeromean_log(zeromean_logSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_estim_lwcc_cpp(data, shrink_int, zeromean_log));
    return rcpp_result_gen;
END_RCPP
}
// sigma_estim_lwnl_cpp
List sigma_estim_lwnl_cpp(arma::mat data, double bandwidth_speed, bool zeromean_log);
RcppExport SEXP _CovEstim_sigma_estim_lwnl_cpp(SEXP dataSEXP, SEXP bandwidth_speedSEXP, SEXP zeromean_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type bandwidth_speed(bandwidth_speedSEXP);
    Rcpp::traits::input_parameter< bool >::type zeromean_log(zeromean_logSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_estim_lwnl_cpp(data, bandwidth_speed, zeromean_log));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CovEstim_sigma_estim_lwident_cpp", (DL_FUNC) &_CovEstim_sigma_estim_lwident_cpp, 3},
    {"_CovEstim_sigma_estim_lwone_cpp", (DL_FUNC) &_CovEstim_sigma_estim_lwone_cpp, 3},
    {"_CovEstim_sigma_estim_lwcc_cpp", (DL_FUNC) &_CovEstim_sigma_estim_lwcc_cpp, 3},
    {"_CovEstim_sigma_estim_lwnl_cpp", (DL_FUNC) &_CovEstim_sigma_estim_lwnl_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_CovEstim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
