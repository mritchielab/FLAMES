#ifndef FLEXIPLEX_H
#define FLEXIPLEX_H

#include <Rcpp.h>

// [[Rcpp::export]]
int flexiplex_cpp(Rcpp::String reads_in, Rcpp::String barcodes_file,
    bool bc_as_readid, int max_bc_editdistance,
    int max_flank_editdistance, Rcpp::StringVector pattern,
    Rcpp::String reads_out, Rcpp::String stats_out,
    Rcpp::String bc_out, int n_threads);

#endif // FLEXIPLEX_H