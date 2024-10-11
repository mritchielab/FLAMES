#ifndef FLEXIPLEX_H
#define FLEXIPLEX_H

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerVector flexiplex_cpp(Rcpp::StringVector reads_in, Rcpp::String barcodes_file,
    bool bc_as_readid, int max_bc_editdistance,
    int max_flank_editdistance, Rcpp::StringVector pattern,
    Rcpp::String reads_out, Rcpp::String stats_out,
    Rcpp::String bc_out, bool reverseComplement, int n_threads);

#endif // FLEXIPLEX_H
