#' Pipeline for Bulk Data
#'
#' Semi-supervised isofrom detection and annotation from long read data.
#' This variant is meant for bulk samples. Specific parameters relating to
#' analysis can be changed either through function arguments, or through a
#' configuration JSON file.
#'
#' @details The default parameters can be changed either through the function
#' arguments are through the configuration JSON file `config_file`. the `pipeline_parameters`
#' section specifies which steps are to be executed in the pipeline - by default, all
#' steps are executed. The `isoform_parameters` section affects isoform detection - key
#' parameters include:
#' * `Min_sup_cnt`, which causes transcripts with less reads aligned than
#' it's value to be discarded;
#' * `MAX_TS_DIST` which merges transcripts with the same intron
#' chain and TSS/TES distace less than `MAX_TS_DIST`;
#' * `strand_specific` which specifies if
#' reads are in the same strand as the mRNA (1), or the reverse complemented (-1) or not
#' strand specific (0), which results in strand information being based on reference annotation.
#'
#' @param annot gene annotations file in gff3  format
#'
#' @param fastq_dir the directory containing fastq files, each containing data from one sample
#'
#' @param in_bam aligned bam file (sorted and indexed). If supplied, this overwrites
#' `fastq_dir` and skips the first alignment step
#'
#' @param outdir directory to store all output files.
#'
#' @param genome_fa genome fasta file.
#'
#' @param minimap2_dir directory containing minimap2, k8 and paftools.js program.
#' k8 and paftools.js are used to convert gff3 to bed12.
#'
#' @param downsample_ratio downsampling ratio if performing downsampling analysis.
#'
#' @param config_file JSON configuration file. If specified, `config_file` overrides
#' all configuration parameters
#'
#' @param do_genome_align Boolean. Specifies whether to run the genome alignment step. `TRUE` is recommended
#' @param do_isoform_id Boolean. Specifies whether to run the isoform identification step. `TRUE` is recommended
#' @param do_read_realign Boolean. Specifies whether to run the read realignment step. `TRUE` is recommended
#' @param do_transcript_quanti Boolean. Specifies whether to run the transcript quantification step. `TRUE` is recommended
#' @param gen_raw-isoform Boolean.
#' @param has_UMI Boolean. Speficies if each gene as a UMI.
#' @param MAX_DIST Numeric
#' @param MAX_TS_DIST Numeric.
#' @param MAX_SPLICE_MATCH_DIST Numeric.
#' @param min_fl_exon_len Numeric.
#' @param Max_site_per_splice Numeric.
#' @param Min_sup_cnt Numeric.
#' @param Min_cnt_pct Numeric.
#' @param Min_sup_pct Numeric.
#' @param strand_specific. 1, -1 or 0. 1 indicates if reads are in the same
#' strand as mRNA, -1 indicates reads are reverse complemented, 0 indicates
#' reads are not strand specific.
#' @param remove_incomp_reads Numeric.
#' @param use_junctions Boolean.
#' @param no_flank Boolean.
#' @param use_annotation Boolean.
#' @param min_tr_coverage Numeric.
#' @param min_read_coverage Numeric.
#'
#' @export
bulk_long_pipeline <- function(annot, fastq_dir, in_bam=NULL, outdir, genome_fa,
                                minimap2_dir=NULL, downsample_ratio=1, config_file=NULL,
                                do_genome_align=TRUE, do_isoform_id=TRUE,
                                do_read_realign=TRUE, do_transcript_quanti=TRUE,
                                gen_raw_isoform=TRUE, has_UMI=FALSE,
                                MAX_DIST=10, MAX_TS_DIST=100, MAX_SPLICE_MATCH_DIST=10,
                                min_fl_exon_len=40, Max_site_per_splice=3, Min_sup_cnt=10,
                                Min_cnt_pct=0.01, Min_sup_pct=0.2, strand_specific=1, remove_incomp_reads=5,
                                use_junctions=TRUE, no_flank=TRUE,
                                use_annotation=TRUE, min_tr_coverage=0.75, min_read_coverage=0.75) {

    # filenames for internal steps
    infq <- paste(outdir, "merged.fastq.gz", sep="/")
    bc_file <- paste(outdir, "pseudo_barcode_annotation.csv", sep="/")

    #if (is.null(in_bam)) {
    # this preprocessing needs only be done if we are using a fastq_dir, instead
    # of a bam file for reads,
    cat("Preprocessing bulk fastqs...\n")
    # create output directory if one doesn't exist
    if (!dir.exists(outdir)) {
        cat("Output directory does not exists: one is being created\n")
        dir.create(outdir)
        print(outdir)
    }
    # run the merge_bulk_fastq function as preprocessing
    merge_bulk_fastq(fastq_dir, bc_file, infq)
    #} else {
    #    bc_file = NULL;
    #}
    generic_long_pipeline(annot, infq, in_bam, outdir, genome_fa,
                minimap2_dir, downsample_ratio, config_file,
                do_genome_align, do_isoform_id,
                do_read_realign, do_transcript_quanti,
                gen_raw_isoform, has_UMI,
                MAX_DIST, MAX_TS_DIST, MAX_SPLICE_MATCH_DIST,
                min_fl_exon_len, Max_site_per_splice, Min_sup_cnt,
                Min_cnt_pct, Min_sup_pct, strand_specific, remove_incomp_reads,
                use_junctions, no_flank,
                use_annotation, min_tr_coverage, min_read_coverage,
                bc_file);

    ## after all this, we need to read in the output data and return R usable datatypes
}
