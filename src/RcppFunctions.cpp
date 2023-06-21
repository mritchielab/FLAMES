#include <string>
#include <vector>
#include <unordered_map>

#include <Rcpp.h>

//' @export
// [[Rcpp::export]]
int flexiplex(Rcpp::String reads_in, Rcpp::String barcodes_file, bool bc_as_readid, int max_bc_editdistance,
              int max_flank_editdistance, Rcpp::List pattern, Rcpp::String reads_out, Rcpp::String stats_out,
              Rcpp::String bc_out, int n_threads);

//' @import zlibbioc
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
