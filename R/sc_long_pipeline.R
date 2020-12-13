#' Pipeline for Single Cell Data
#'
#' DESC
#' 
#' @inheritParams bulk_long_pipeline
#' @export
sc_long_pipeline <- function(annot, infq, in_bam=NULL, outdir, genome_fa,
                                minimap2_dir=NULL, downsample_ratio=1, config_file=NULL,
                                do_genome_align=TRUE, do_isoform_id=TRUE, 
                                do_read_realign=TRUE, do_transcript_quanti=TRUE,
                                gen_raw_isoform=TRUE, has_UMI=FALSE,
                                MAX_DIST=10, MAX_TS_DIST=100, MAX_SPLICE_MATCH_DIST=10,
                                min_fl_exon_len=40, Max_site_per_splice=3, Min_sup_cnt=10, 
                                Min_cnt_pct=0.01, Min_sup_pct=0.2, strand_specific=1, remove_incomp_reads=5,
                                use_junctions=TRUE, no_flank=TRUE,
                                use_annotation=TRUE, min_tr_coverage=0.75, min_read_coverage=0.75) {
    # this seems to be the same function as bulk_long_pipeline, except without the mergine
    # of fastq files at the start.

        # seems like we need an extra step in here before we do this.
        # something involving a new function match_cell_barcode
    generic_long_pipeline(annot, infq, in_bam, outdir, genome_fa, 
            minimap2_dir, downsample_ratio, config_file,
            do_genome_align, do_isoform_id,
            do_read_realign, do_transcript_quanti,
            gen_raw_isoform, has_UMI,
            MAX_DIST, MAX_TS_DIST, MAX_SPLICE_MATCH_DIST,
            min_fl_exon_len, Max_site_per_splice, Min_sup_cnt,
            Min_cnt_pct, Min_sup_pct, strand_specific, remove_incomp_reads,
            use_junctions, no_flank,
            use_annotation, min_tr_coverage, min_read_coverage);
}

