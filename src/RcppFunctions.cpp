#include <string>
#include <vector>
#include <unordered_map>

#include <Rcpp.h>

#include "main-functions/match_cell_barcode.h"

//' Match Cell Barcodes
//'
//' @description Match cell barcodes in the given fastq directory with the reference csv, \code{ref_csv}. Matches are returned
//' in the output file \code{out_fastq}
//' The flanking sequence is aligned to the first 30000 reads to identify the regions where cell barcode is likely to be found within. 
//' Next, sequences within this region are matched to barcodes in \code{ref_csv}, allowing \code{MAX_DIST} hamming distances. 
//' Reads that are successfully matched with a barcode are reported as the \code{barcode hm match} count. 
//' Every read that could not be matched in the previous step is aligned to the flanking sequence again to identify the location of 
//' barcode individually, and barcode matching is done with up to \code{MAX_DIST} levenshtein distances (allowing indels). Reads that are
//' matched by this step is reported as the \code{fuzzy match} counts.
//'
//' @param fastq_dir directory containing fastq files to match
//' @param stats_file NEEDED
//' @param out_fastq output filename for matched barcodes
//' @param ref_csv NEEDED
//' @param MAX_DIST int; maximum edit distance
//' @param UMI_LEN int; length of UMI sequences
//'
//' @return returns NULL
//' @import zlibbioc
//' @useDynLib FLAMES, .registration=TRUE
//' @export
// [[Rcpp::export]]
void
match_cell_barcode_cpp
(
    Rcpp::String fastq_dir, 
    Rcpp::String stats_file, 
    Rcpp::String out_fastq, 
    Rcpp::String ref_csv, 
    int MAX_DIST, 
    int UMI_LEN = 10
)
{
    return match_cell_barcode(
        fastq_dir,
        stats_file,
        out_fastq,
        ref_csv,
        MAX_DIST,
        UMI_LEN
    );
}

/*****************************************************************************/

#include "main-functions/merge_bulk.h"
//' @useDynLib FLAMES, .registration=TRUE
// [[Rcpp::export]]
void 
merge_bulk_fastq_cpp(Rcpp::StringVector fastq_files, Rcpp::String out_fastq) 
{
    return merge_bulk_fastq(fastq_files, out_fastq);
}
