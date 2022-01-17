#include <Rcpp.h>

// a file for holding a bunch of wrapper functions
// until we figure out how to get Rcpp looking in subdirectories

/*****************************************************************************/

#include "file-handling/parse_json_config.h"

// [[Rcpp::export]]
Rcpp::List
parse_json_config_cpp (std::string jsonFile)
{
    return parse_json_config(jsonFile);
}

// [[Rcpp::export]]
void
print_config_cpp (Rcpp::List config)
{
    return print_config(config);
}

/*****************************************************************************/

#include "file-handling/gtf_to_bed.h"

// [[Rcpp::export]]
void
gtf_to_bed_cpp
(
    std::string in_gtf, 
    std::string out_bed,
    std::string chrom_sizes_file
)
{
    return gtf_to_bed(in_gtf, out_bed, chrom_sizes_file);
}

/*****************************************************************************/

#include "file-handling/parse_realigned_bam.h"

// [[Rcpp::export]]
void
read_entire_bam_cpp
(
    std::string bam_in,
    std::string log_out
)
{
    return read_entire_bam(bam_in, log_out);
}

/*****************************************************************************/

#include "main-functions/find_isoform.h"

// [[Rcpp::export]]
Rcpp::List
find_isoform_cpp
(
    std::string gff3,
    std::string genome_bam,
    std::string isoform_gff3,
    std::string tss_test_stat,
    std::string genomefa,
    std::string transcript_fa,
    int         downsample_ratio,
    Rcpp::List  config_list,
    std::string raw_splice_isoform
)
{
    return find_isoform(
        gff3, 
        genome_bam,
        isoform_gff3,
        tss_test_stat,
        genomefa,
        transcript_fa,
        downsample_ratio,
        config_list,
        raw_splice_isoform);
}

/*****************************************************************************/

#include "main-functions/group_bam2isoform.h"

// [[Rcpp::export]]
void
bam_read_cpp
(
    std::string bam_in, 
    std::string chr, 
    int s, 
    int e
)
{
    return bam_read(bam_in, chr, s, e);
}

/*****************************************************************************/

#include "main-functions/match_cell_barcode.h"

// [[Rcpp::export]]
void
match_cell_barcode_cpp
(
    String fastq_dir, 
    String stats_file, 
    String out_fastq, 
    String ref_csv, 
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
merge_bulk_fastq_cpp(StringVector fastq_files, String out_fastq) 
{
    return merge_bulk_fastq(fastq_files, out_fastq);
}

/*****************************************************************************/

#include "main-functions/quantification.h"

// [[Rcpp::export]]
void
quantification_cpp
(
    Rcpp::List  config_list,
    std::string realign_bam,
    std::string transcript_fa_idx,
    Rcpp::List  isoform_objects_list,
    std::string tr_cnt_csv,
    std::string tr_badcov_cnt_csv,
    std::string isoform_gff3,
    std::string annot,
    std::string isoform_gff3_f,
    std::string FSM_anno_out
)
{
    quantification(
        config_list,
        realign_bam,
        transcript_fa_idx,
        isoform_objects_list,
        tr_cnt_csv,
        tr_badcov_cnt_csv,
        isoform_gff3,
        annot,
        isoform_gff3_f,
        FSM_anno_out
    );
}

/*****************************************************************************/

#include "utility/minimap2/minimap2_align.h"

// [[Rcpp::export]]
void
minimap2_align_cpp
(
    std::string mm2_prog_path,
    std::string fa_file,
    std::string fq_in,
    std::string sam_out,
    bool no_flank,
    std::string bed12_junc
)
{
    return minimap2_align(
        mm2_prog_path,
        fa_file,
        fq_in,
        sam_out,
        no_flank,
        bed12_junc);
}