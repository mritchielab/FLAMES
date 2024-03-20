#ifndef PILEUP_H
#define PILEUP_H

// [[Rcpp::export]]
Rcpp::NumericMatrix
variant_count_matrix_cpp(Rcpp::String bam_path, Rcpp::String seqname, int pos,
                         bool indel, Rcpp::StringVector barcodes, bool verbose);

#endif // PILEUP_H
