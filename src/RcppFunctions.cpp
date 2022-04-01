#include <string>
#include <vector>
#include <unordered_map>

#include <Rcpp.h>


#include "main-functions/match_cell_barcode.h"

//' Match Cell Barcodes
//'
//' @description Match cell barcodes in the given fastq directory with the reference csv, \code{ref_csv}. Matches are returned
//' in the output file \code{out_fastq}
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

// [[Rcpp::export]]
void 
merge_bulk_fastq_cpp(Rcpp::StringVector fastq_files, Rcpp::String out_fastq) 
{
    return merge_bulk_fastq(fastq_files, out_fastq);
}


/*****************************************************************************/
// GffRead
// https://github.com/gpertea/gffread MIT License
// https://github.com/gpertea/gclib Artistic License 2.0
// TODO: fix anonymous structs (set options in Makevar?)

#include "gffread/gffread.h"

// [[Rcpp::export]]
int
gffread_cpp(
    Rcpp::String genome_fa,
    Rcpp::String transcript_fa,
    Rcpp::String gff3
)
{
    int argc = 6;
    char* argv[] = {(char *)"gffread", (char*)gff3.get_cstring(), (char *)"-g", (char*)genome_fa.get_cstring(), (char *)"-w", (char*)transcript_fa.get_cstring()};
    return gffread(argc, argv);
}