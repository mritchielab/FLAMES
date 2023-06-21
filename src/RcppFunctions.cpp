#include <string>
#include <vector>
#include <unordered_map>

#include <Rcpp.h>

//' @export
// [[Rcpp::export]]
int flexiplex(Rcpp::String reads_in, Rcpp::String barcodes_file, bool bc_as_readid, int max_bc_editdistance,
              int max_flank_editdistance, Rcpp::List pattern, Rcpp::String reads_out, Rcpp::String stats_out,
              Rcpp::String bc_out, int n_threads);
