// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// flames_test_func
int flames_test_func();
RcppExport SEXP _FLAMES_flames_test_func() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(flames_test_func());
    return rcpp_result_gen;
END_RCPP
}
// bam_read
void bam_read(std::string bam_in);
RcppExport SEXP _FLAMES_bam_read(SEXP bam_inSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bam_in(bam_inSEXP);
    bam_read(bam_in);
    return R_NilValue;
END_RCPP
}
// match_cell_barcode
void match_cell_barcode(String fastq_dir, String stats_file, String out_fastq, String ref_csv, int MAX_DIST, int UMI_LEN);
RcppExport SEXP _FLAMES_match_cell_barcode(SEXP fastq_dirSEXP, SEXP stats_fileSEXP, SEXP out_fastqSEXP, SEXP ref_csvSEXP, SEXP MAX_DISTSEXP, SEXP UMI_LENSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type fastq_dir(fastq_dirSEXP);
    Rcpp::traits::input_parameter< String >::type stats_file(stats_fileSEXP);
    Rcpp::traits::input_parameter< String >::type out_fastq(out_fastqSEXP);
    Rcpp::traits::input_parameter< String >::type ref_csv(ref_csvSEXP);
    Rcpp::traits::input_parameter< int >::type MAX_DIST(MAX_DISTSEXP);
    Rcpp::traits::input_parameter< int >::type UMI_LEN(UMI_LENSEXP);
    match_cell_barcode(fastq_dir, stats_file, out_fastq, ref_csv, MAX_DIST, UMI_LEN);
    return R_NilValue;
END_RCPP
}
// merge_bulk_fastq_cpp
void merge_bulk_fastq_cpp(StringVector fastq_files, String out_fastq);
RcppExport SEXP _FLAMES_merge_bulk_fastq_cpp(SEXP fastq_filesSEXP, SEXP out_fastqSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type fastq_files(fastq_filesSEXP);
    Rcpp::traits::input_parameter< String >::type out_fastq(out_fastqSEXP);
    merge_bulk_fastq_cpp(fastq_files, out_fastq);
    return R_NilValue;
END_RCPP
}
// parse_json_config_cpp
Config parse_json_config_cpp(std::string json_file);
RcppExport SEXP _FLAMES_parse_json_config_cpp(SEXP json_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type json_file(json_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(parse_json_config_cpp(json_file));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FLAMES_flames_test_func", (DL_FUNC) &_FLAMES_flames_test_func, 0},
    {"_FLAMES_bam_read", (DL_FUNC) &_FLAMES_bam_read, 1},
    {"_FLAMES_match_cell_barcode", (DL_FUNC) &_FLAMES_match_cell_barcode, 6},
    {"_FLAMES_merge_bulk_fastq_cpp", (DL_FUNC) &_FLAMES_merge_bulk_fastq_cpp, 2},
    {"_FLAMES_parse_json_config_cpp", (DL_FUNC) &_FLAMES_parse_json_config_cpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_FLAMES(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
