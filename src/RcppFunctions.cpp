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
//' @examples
//' outdir <- tempfile()
//' dir.create(outdir)
//' bc_allow <- file.path(outdir, "bc_allow.tsv")
//' R.utils::gunzip(filename = system.file("extdata/bc_allow.tsv.gz", package = "FLAMES"), destname = bc_allow, remove = FALSE)
//' find_barcode(
//'    fastq_dir = system.file("extdata/fastq", package = "FLAMES"),
//'    stats_file = file.path(outdir, "bc_stat"),
//'    out_fastq = file.path(outdir, "demultiplexed.fq.gz"),
//'    ref_csv = bc_allow,
//'    MAX_DIST = 2,
//'    UMI_LEN = 10
//')
//' @export
// [[Rcpp::export]]
void
find_barcode
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

#include "main-functions/find_isoform.h"

// [[Rcpp::export]]
void
find_isoform_multithread
(
    const std::string &gff3,
    const std::string &genome_bam,
    const std::string &isoform_gff3,
    const std::string &tss_tes_stat,
    const std::string &genomefa,
    const std::string &transcript_fa,
    const Rcpp::List  &isoform_parameters,
    const std::string &raw_splice_isoform
)
{
    find_isoform_multithread_cpp(
        gff3, 
        genome_bam,
        isoform_gff3,
        tss_tes_stat,
        genomefa,
        transcript_fa,
        isoform_parameters,
        raw_splice_isoform
    );
}
