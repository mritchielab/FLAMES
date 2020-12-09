#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @importFrom reticulate import_from_path
#' @export
sc_long_pipeline <- function(annot, infq, inbam=NULL, outdir, genomefa, 
                    minimap2_dir=NULL, config_file=NULL, downsample_ratio=1) {

    if (is.null(config_file)) { #check these are the right parameters
        config = 
            list(
                pipeline_parameters=
                    list(do_genome_alignment=do_genome_align, do_isoform_identification=do_isoform_id, 
                        do_read_realignment=do_read_realign, do_transcript_quantification=do_transcript_quantification),
                global_parameters=
                    list(generate_raw_isoform=gen_raw_isoform, has_UMI=has_UMI),
                isoform_parameters=
                    list(MAX_DIST=MAX_DIST, MAX_TS_DIST=MAX_TS_DIST, MAX_SPLICE_MATCH_DIST=MAX_SPLICE_MATCH_DIST,
                                min_fl_exon_len=min_fl_exon_len, Max_site_per_splice=Max_site_per_splice, Min_sup_cnt=Min_sup_cnt, 
                                Min_cnt_pct=Min_cnt_pct, Min_sup_pct=Min_sup_pct.2, strand_specific=strand_specific, 
                                remove_incomp_reads=remove_incomp_reads),
                alignment_parameters=
                    list(use_junctions=use_junctions, no_flank=no_flank),
                realign_parameters=
                    list(use_annotation=use_annotation),
                transcript_counting=
                    list(min_tr_coverage=min_tr_coverage, min_read_coverage=min_read_coverage) 
                )
        if (MAX_DIST <= 0 || MAX_TS_DIST <= 0 || MAX_SPLICE_MATCH_DIST <= 0 || Max_site_per_splce <= 0 ||
            Min_sup_cnt <= 0 || Min_sup_pct <= 0 || (strand_specific != -1 && strand_specific != 0 && strand_specific != 1) || remove_incomp_reads < 0) {
                stop("MAX_DIST,  MAX_TS_DIST, MAX_SPLICE_MATCH_DIST,  Max_site_per_splce, Min_sup_cnt and Min_sup_pct must be greater than 0. strand_specific must be -1, 0 or 1 and remove_incomp_reads must be >= 0.")
        }
    } else {
        config = parse_json_config(config_file)
    }

    # check non config argument
    if (downsample_ratio > 1 || downsample_ratio <= 0) {
        stop("downsample_ratio should be between 0 and 1")
    }
    if (!file.exists(infq) || !file.exists(annot) || !file.exists(genome_fa)) {
        stop(paste0("Make sure all files exists: ", infq, ", ", annot, ", ", genome_fa))
    }
    if (is.null(in_bam)) {
        in_bam = ""
    } else {
        if (!file.exists(in_bam)) {
            stop("Make sure input in_bam file exists")
        }
    }
    if (is.null(minimap2_dir)) minimap2_dir = ""

    isoform_gff3 = paste(outdir, "isoform_annotated.gff3", sep="/")
    isoform_gff3_f = paste(outdir, "isoform_annotated.filtered.gff3", sep="/")
    FSM_anno_out = paste(outdir, "isoform_FSM_annotation.csv", sep="/")
    raw_splice_isoform = paste(outdir, "splice_raw.gff3", sep="/")
    tss_tes_stat = paste(outdir, "tss_tes.bedgraph", sep="/")
    transcript_fa = paste(outdir, "transcript_assembly.fa", sep="/")
    transcript_fa_idx = paste(outdir, "transcript_assembly.fa.fai", sep="/")
    tmp_bam = paste(outdir, "tmp.align.bam", sep="/")
    tmp_bed = paste(outdir, "tmp.splice_anno.bed12", sep="/")
    genome_bam = paste(outdir, "align2genome.bam", sep="/")
    realign_bam = paste(outdir, "realign2transcript.bam", sep="/")
    tr_cnt_csv = paste(outdir, "transcript_count.csv.gz", sep="/")
    tr_badcov_cnt_csv = paste(outdir, "transcript_count.bad_coverage.csv.gz", sep="/")

}

