// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// adapt
List adapt(const arma::mat& mydata, int iter_times);
RcppExport SEXP _RARE_adapt(SEXP mydataSEXP, SEXP iter_timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mydata(mydataSEXP);
    Rcpp::traits::input_parameter< int >::type iter_times(iter_timesSEXP);
    rcpp_result_gen = Rcpp::wrap(adapt(mydata, iter_times));
    return rcpp_result_gen;
END_RCPP
}
// calculateP
double calculateP(double mean, double sd);
RcppExport SEXP _RARE_calculateP(SEXP meanSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(calculateP(mean, sd));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _RARE_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RARE_adapt", (DL_FUNC) &_RARE_adapt, 2},
    {"_RARE_calculateP", (DL_FUNC) &_RARE_calculateP, 2},
    {"_RARE_rcpp_hello_world", (DL_FUNC) &_RARE_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_RARE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
