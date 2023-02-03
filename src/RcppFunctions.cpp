#include <string>
#include <vector>
#include <unordered_map>

#include <Rcpp.h>

#include "main-functions/match_cell_barcode.h"

//' @import zlibbioc
//' @useDynLib FLAMES, .registration=TRUE
// [[Rcpp::export]]
Rcpp::List
match_cell_barcode
(
    Rcpp::String fastq_dir, 
    Rcpp::String stats_file, 
    Rcpp::String out_fastq, 
    Rcpp::String ref_csv, 
    int MAX_DIST, 
    int UMI_LEN,
    Rcpp::String left_seq,
    int min_length,
    bool reverse_complement,
    bool fixed_range
);
