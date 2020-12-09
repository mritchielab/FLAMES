#' Title
#'
#' Semi-supervised isofrom detection and annotation from long read data.
#' This variant is meant for bulk samples. THIS FUNCTION IS GOING TO BE THE ENTIRE 
#' FLAMES PROCESS. IT SHOULD BE BUILT FROM THE SMALLER FUNCTIONS THAT PREVIOUSLY WERE
#' CARRIED OUT UNDER THE HOOD. MOST FUNCTIONS THAT THIS ONE USES SHOULD BE AVAILABLE TO THE USER.
#'
#'  THIS NEED PROPER TYPES
#' @param annot gene annotations in gff3 file format. str
#'
#' @param fastq_dir the folder containing fastq files, each containing data from one sample. str
#'
#' @param in_bam aligned bam file (sorted and indexed). Overwrite --infq? str
#'
#' @param outdir directory to deposite all results in rootdir. Use absolute path? str
#'
#' @param genome_fa genome fasta file. str
#'
#' @param minimap2_dir directory containing minimap2, k8 and paftools.js program. k8 and paftools.js are used to convert gff3 to bed12. str
#'
#' @param config_file json configuration files. str
#'
#' @param downsample_ratio downsampling ratio if performing downsampling analysis. str
#' @importFrom reticulate import_from_path
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
            
            # should we just make them all optional parameters with default values, and give
            # and optional config JSON file which will override the defaults?
    # setup config file if none is given, or read in json config
    if (is.null(config_file)) {
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

    # filenames for internal steps
    infq <- paste(outdir, "merged.fastq.gz", sep="/")
    bc_file <- paste(outdir, "pseudo_barcode_annotation.csv", sep="/")

    cat("Preprocessing bulk fastqs...\n")
    # create output directory if one doesn't exist
    if (!dir.exists(outdir)) {
        cat("Output directory does not exists: one is being created\n")
        dir.create(outdir)
        print(outdir)
    }
    # run the merge_bulk_fastq function as preprocessing
    merge_bulk_fastq(fastq_dir, bc_file, infq)

    cat("Running FLAMES pipeline...\n")
    # argument verificiation
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

    # setup of internal arguments which hold output files and intermediate files
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

    cat("#### Input parameters:\n")
    print_config(config)
    cat("\tgene annotation:", annot, "\n")
    cat("\tgenome fasta:", genome_fa, "\n")
    if (in_bam != "") {
        cat("\tinput bam:", in_bam, "\n")
        genome_bam = in_bam
    } else cat("\tinput fastq:", infq, "\n")
    cat("\toutput directory:", outdir, "\n")
    cat("\tdirectory containing minimap2:", minimap2_dir, "\n")

    # align reads to genome
    if (in_bam=="" && config$pipeline_parameters$do_genome_alignment) {
        cat("#### Aligning reads to genome using minimap2\n")

        tmp_bed <- paste(outdir, "tmp.splice_anno.bed12", sep=.Platform$file.sep)
        tmp_bam <- paste(outdir, "tmp.align.bam", sep=.Platform$file.sep)

        if (config$alignment_parameters$use_junctions) {
            gff3_to_bed12(minimap2_dir, annot, tmp_bed)
        }
        minimap2_align(minimap2_dir, genome_fa, infq, tmp_bam, 
            no_flank=config$alignment_parameters$no_flank, bed12_junc = if (config$alignment_parameters$use_junctions) tmp_bed else NULL)
        samtools_sort_index(tmp_bam, genome_bam)
        file.remove(tmp_bam)
        if (config$alignment_parameters$use_junctions) file.remove(tmp_bed)
    } else {
        cat("#### Skip aligning reads to genome\n")
    }

    # find isofrom
    cat("#### Read genne annotations\n")
    gff3_parse_result <- parse_gff_tree(annot)
    chr_to_gene = gff3_parse_result$chr_to_gene # chr_to_gene appears to be the exact same as python returns, so not sure what the issue is
    transcript_dict = gff3_parse_result$transcript_dict
    gene_to_transcript = gff3_parse_result$gene_to_transcript
    transcript_to_exon = gff3_parse_result$transcript_to_exon
    remove_similar_tr(transcript_dict, gene_to_transcript, transcript_to_exon)

    if (config$pipeline_parameters$do_isoform_identification) {
        cat("#### Find isoforms\n")
        transcript_to_junctions = list()
        for (tr in names(transcript_to_exon)) {
            transcript_to_junctions[[tr]] = blocks_to_junctions(transcript_to_exon[[tr]])
        }
        gene_dict <- get_gene_flat(gene_to_transcript, transcript_to_exon)
        chr_to_blocks <- get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
        group_bam2isoform(genome_bam, isoform_gff3, tss_tes_stat, "", chr_to_blocks, 
                gene_dict, transcript_to_junctions, transcript_dict, genomefa,
                config=config$isoform_parameters, downsample_ratio=downsample_ratio, 
                raw_gff3=if (config$global_parameters$generate_raw_isoform) raw_splice_isoform else NULL)
    } else {
        ## skip finding isoform.
        cat("#### Skip finding isoforms\n")
    }
    # get fasta
    isoform_gff3_parse <- parse_gff_tree(isoform_gff3)
    chr_to_gene_i <- isoform_gff3_parse$chr_to_gene
    transcript_dict_i <- isoform_gff3_parse$transcript_dict
    gene_to_transcript_i <- isoform_gff3_parse$gene_to_transcript
    transcript_to_exon_i <- isoform_gff3_parse$transcript_to_exon

    if (!config$realign_parameters$use_annotation) gff3_parse_result = NULL
    get_transcript_seq(genome_fa, transcript_fa, chr_to_gene_i, transcript_dict_i, 
            gene_to_transcript_i, transcript_to_exon_i, ref_dict=gff3_parse_result)
    
    # realign to transcript
    if (do_read_realign) {
        cat("#### Realign to transcript using minimap2\n")
        tmp_bam <- paste(outdir, "tmp.align.bam", sep=.Platform$file.sep)
        minimap2_tr_align(minimap2_dir, transcript_fa, infq, tmp_bam)
        samtools_sort_index(tmp_bam, realign_bam)
        file.remove(tmp_bam)
    } else {
        cat("#### Skip read realignment\n")
    }

    #quantification
    if (config$pipeline_parameters$do_transcript_quantification) {
        cat("#### Generating transcript count matrix\n")
        parse_realign <- parse_realigned_bam(realign_bam, transcript_fa_idx,
            config$isoform_parameters$Min_sup_cnt,
            config$transcript_counting$min_tr_coverage,
            config$transcript_counting$min_read_coverage, bc_file=bc_file)
        
        tr_cnt = wrt_tr_to_csv(parse_realign$bc_tr_count_dict, transcript_dict_i, tr_cnt_csv, transcript_dict, config$global_parameters$has_UMI)
        wrt_tr_to_csv(parse_realign$bc_tr_badcov_count_dict, transcript_dict_i, tr_badcov_cnt_csv, transcript_dict, config$global_parameters$has_UMI)
        annotate_filter_gff(isoform_gff3, annot, isoform_gff3_f, FSM_anno_out, tr_cnt, config$isoform_parameters$Min_sup_cnt)
    } else {
        cat("#### Skip transcript quantification\n")
    }

    ## after all this, we need to read in the output data and return R usable datatypes
}