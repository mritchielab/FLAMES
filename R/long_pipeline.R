#' Generic FLAMES pipeline
#'
#' Generic implementation of the flames pipeline. Used for both bulk reads and
#' single cell reads.
#'
#' @inheritParams sc_long_pipeline
#' @param in_bam optional BAM file which replaces fastq directory argument. This skips the genome alignment and
#' realignment steps
#' @param bc_file file containing the pseudo barcode annotations generated
#' from bulk_long_pipeline. If given, it is used for quantification.
#' 
#' @return This generic function returns NULL, instead providing output files
#' in the given `outdir` directory. These files are loaded into R in either
#' a SummarizedExperiment or SingleCellExperiment object by the callers to this
#' function, `sc_long_pipeline()` and `bulk_long_pipeline()` respectively.
generic_long_pipeline <- function(annot, fastq, in_bam, outdir, genome_fa,
                minimap2_dir, downsample_ratio, config_file,
                do_genome_align, do_isoform_id=TRUE,
                do_read_realign, do_transcript_quanti,
                gen_raw_isoform, has_UMI,
                MAX_DIST, MAX_TS_DIST, MAX_SPLICE_MATCH_DIST,
                min_fl_exon_len, Max_site_per_splice, Min_sup_cnt,
                Min_cnt_pct, Min_sup_pct, strand_specific, remove_incomp_reads,
                use_junctions, no_flank,
                use_annotation, min_tr_coverage, min_read_coverage,
                bc_file=NULL) {

    if (!dir.exists(outdir)) {
        cat("Output directory does not exists: one is being created\n")
        dir.create(outdir)
        print(outdir)
    }
    if (is.null(config_file)) {
        config = create_config(do_genome_align, do_isoform_id,
                        do_read_realign, do_transcript_quanti,
                        gen_raw_isoform, has_UMI,
                        MAX_DIST, MAX_TS_DIST, MAX_SPLICE_MATCH_DIST,
                        min_fl_exon_len, Max_site_per_splice, Min_sup_cnt,
                        Min_cnt_pct, Min_sup_pct, strand_specific, remove_incomp_reads,
                        use_junctions, no_flank,
                        use_annotation, min_tr_coverage, min_read_coverage)
    } else {
        config = parse_json_config(config_file)
    }

    if (!config$pipeline_parameters$do_isoform_identification) {
        stop("Isoform Identification is required for FLAMES execution. Change this value in the configuration file or \
        set the argument as TRUE. If isoform identification is not required, you can manually execute the pipeline by following \
        the vignette")
    }

    cat("Running FLAMES pipeline...\n")
    # argument verificiation
    if (downsample_ratio > 1 || downsample_ratio <= 0) {
        stop("downsample_ratio should be between 0 and 1")
    }
    if (!is.null(fastq) && !file.exists(fastq)) stop(paste0("Make sure ", fastq, " exists."))
    if (!file.exists(annot)) stop(paste0("Make sure ", annot, " exists."))
    if (!file.exists(genome_fa)) stop(paste0("Make sure ", genome_fa, " exists."))

    using_bam = FALSE
    if (!is.null(in_bam)) {
        using_bam = TRUE
        bc_file = NULL
        fastq = NULL
        if (!file.exists(in_bam)) stop ("Make sure in_bam exists")
    } 
    if (is.null(minimap2_dir)) minimap2_dir = ""

    # setup of internal arguments which hold output files and intermediate files
    isoform_gff3 = paste(outdir, "isoform_annotated.gff3", sep="/")
    isoform_gff3_f = paste(outdir, "isoform_annotated.filtered.gff3", sep="/")
    FSM_anno_out = paste(outdir, "isoform_FSM_annotation.csv", sep="/")
    raw_splice_isoform = paste(outdir, "splice_raw.gff3", sep="/")
    tss_tes_stat = paste(outdir, "tss_tes.bedgraph", sep="/")
    transcript_fa = paste(outdir, "transcript_assembly.fa", sep="/")
    transcript_fa_idx = paste(outdir, "transcript_assembly.fa.fai", sep="/")
    tmp_bam = paste(outdir, "tmp_align.bam", sep="/")
    tmp_bed = paste(outdir, "tmp_splice_anno.bed12", sep="/")
    tmp_sam = paste(outdir, "tmp_align.sam", sep="/")
    genome_bam = paste(outdir, "align2genome.bam", sep="/")
    realign_bam = paste(outdir, "realign2transcript.bam", sep="/")
    tr_cnt_csv = paste(outdir, "transcript_count.csv.gz", sep="/")
    tr_badcov_cnt_csv = paste(outdir, "transcript_count.bad_coverage.csv.gz", sep="/")

    cat("#### Input parameters:\n")
    print_config(config)
    cat("\tgene annotation:", annot, "\n")
    cat("\tgenome fasta:", genome_fa, "\n")
    if (using_bam) {
        cat("\tinput bam:", in_bam, "\n")
        genome_bam = in_bam
    } else cat("\tinput fastq:", fastq, "\n")
    cat("\toutput directory:", outdir, "\n")
    cat("\tdirectory containing minimap2:", minimap2_dir, "\n")

    # align reads to genome
    #if (!using_bam && config$pipeline_parameters$do_genome_alignment) {
    if (config$pipeline_parameters$do_genome_alignment) {
        cat("#### Aligning reads to genome using minimap2\n")

        if (config$alignment_parameters$use_junctions) {
            gff3_to_bed12(minimap2_dir, annot, tmp_bed)
        }
        minimap2_align(minimap2_dir, genome_fa, fastq, tmp_sam,
            no_flank=config$alignment_parameters$no_flank, bed12_junc = if (config$alignment_parameters$use_junctions) tmp_bed else NULL)
        samtools_as_bam(tmp_sam, tmp_bam)
        samtools_sort_index(tmp_bam, genome_bam)
        file.remove(tmp_sam)
        file.remove(tmp_bam)
        if (config$alignment_parameters$use_junctions) file.remove(tmp_bed)
    } else {
        cat("#### Skip aligning reads to genome\n")
    }

    # find isofrom
    cat("#### Read genne annotations\n")
    gff3_parse_result <- parse_gff_tree(annot)
    chr_to_gene = gff3_parse_result$chr_to_gene
    transcript_dict = gff3_parse_result$transcript_dict
    gene_to_transcript = gff3_parse_result$gene_to_transcript
    transcript_to_exon = gff3_parse_result$transcript_to_exon
    remove_similar_tr(gene_to_transcript, transcript_to_exon) # issue

    # do_isoform_identification is now always required to be true. 
    #if (config$pipeline_parameters$do_isoform_identification) {
    cat("#### Find isoforms\n")
    transcript_to_junctions = list()
    for (tr in names(transcript_to_exon)) {
        transcript_to_junctions[[tr]] = blocks_to_junctions(transcript_to_exon[[tr]])
    }
    gene_dict <- get_gene_flat(gene_to_transcript, transcript_to_exon)
    chr_to_blocks <- get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    group_bam2isoform(genome_bam, isoform_gff3, tss_tes_stat, "", chr_to_blocks,
            gene_dict, transcript_to_junctions, transcript_dict, genome_fa,
            config=config$isoform_parameters, downsample_ratio=downsample_ratio,
            raw_gff3=if (config$global_parameters$generate_raw_isoform) raw_splice_isoform else NULL)
    #} else {
    #    ## skip finding isoform.
    #    cat("#### Skip finding isoforms\n")
    #}

    # get fasta
    isoform_gff3_parse <- parse_gff_tree(isoform_gff3)
    chr_to_gene_i <- isoform_gff3_parse$chr_to_gene
    transcript_dict_i <- isoform_gff3_parse$transcript_dict ## this is the only variable required after get_transcript_seq
    gene_to_transcript_i <- isoform_gff3_parse$gene_to_transcript
    transcript_to_exon_i <- isoform_gff3_parse$transcript_to_exon

    get_transcript_seq(genome_fa, transcript_fa, chr_to_gene_i, transcript_dict_i,
            gene_to_transcript_i, transcript_to_exon_i, ref_dict=if (config$realign_parameters$use_annotation) gff3_parse_result else NULL)

    # realign to transcript
    #if (!using_bam && do_read_realign) {
    if (do_read_realign) {
        cat("#### Realign to transcript using minimap2\n")
        minimap2_tr_align(minimap2_dir, transcript_fa, fastq, tmp_sam)
        samtools_as_bam(tmp_sam, tmp_bam)
        samtools_sort_index(tmp_bam, realign_bam)
        file.remove(tmp_sam)
        file.remove(tmp_bam)
    } else {
        cat("#### Skip read realignment\n")
    }

    #quantification
    if (config$pipeline_parameters$do_transcript_quantification) {
        cat("#### Generating transcript count matrix\n")
        #if (using_bam) {
            # i have no idea what to do if input is a bam file. This does not work as realign_bam and transcript_fa_idx are required, 
            # which are both produced from minimap2_tr_align.
            #realign_bam <- in_bam
        #}
        if (is.null(bc_file)) {
            # sc_long_pipeline version of realigned bam
            parse_realign <- parse_realigned_bam(realign_bam, transcript_fa_idx,
                config$isoform_parameters$Min_sup_cnt,
                config$transcript_counting$min_tr_coverage,
                config$transcript_counting$min_read_coverage)
        } else {
            parse_realign <- parse_realigned_bam(realign_bam, transcript_fa_idx,
                config$isoform_parameters$Min_sup_cnt,
                config$transcript_counting$min_tr_coverage,
                config$transcript_counting$min_read_coverage,
                bc_file=bc_file)
        }
        tr_cnt = wrt_tr_to_csv(parse_realign$bc_tr_count_dict, transcript_dict_i, tr_cnt_csv, transcript_dict, config$global_parameters$has_UMI)
        wrt_tr_to_csv(parse_realign$bc_tr_badcov_count_dict, transcript_dict_i, tr_badcov_cnt_csv, transcript_dict, config$global_parameters$has_UMI)
        annotate_filter_gff(isoform_gff3, annot, isoform_gff3_f, FSM_anno_out, tr_cnt, config$isoform_parameters$Min_sup_cnt)
    } else {
        cat("#### Skip transcript quantification\n")
    }

}

