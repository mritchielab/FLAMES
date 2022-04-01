#' Generic FLAMES pipeline
#'
#' Generic implementation of the flames pipeline. Used for both bulk reads and
#' single cell reads.
#'
#' @inheritParams sc_long_pipeline
#' @param in_bam optional BAM file which replaces fastq directory argument. This skips the genome alignment and
#' realignment steps
#'
#' @return This generic function returns a named list containing the output file names of the provided output files
#' in the given `outdir` directory. These files are loaded into R in either
#' a SummarizedExperiment or SingleCellExperiment object by the callers to this
#' function, `sc_long_pipeline()` and `bulk_long_pipeline()` respectively.
#'
#' @importFrom GenomeInfoDb seqlengths
generic_long_pipeline <-
    function(annot,
             fastq,
             in_bam,
             outdir,
             genome_fa,
             minimap2_dir,
             seed = NULL,
             downsample_ratio,
             config_file,
             do_genome_align,
             do_isoform_id = TRUE,
             isoform_id_bambu = FALSE,
             do_read_realign,
             do_transcript_quanti,
             gen_raw_isoform,
             has_UMI,
             MAX_DIST,
             MAX_TS_DIST,
             MAX_SPLICE_MATCH_DIST,
             min_fl_exon_len,
             Max_site_per_splice,
             Min_sup_cnt,
             Min_cnt_pct,
             Min_sup_pct,
             strand_specific,
             remove_incomp_reads,
             use_junctions,
             no_flank,
             use_annotation,
             min_tr_coverage,
             min_read_coverage) {
        cat("Running FLAMES pipeline...\n")
        config <- parse_json_config(config_file)

        using_bam <- FALSE
        if (!is.null(in_bam)) {
            if (!file.exists(paste0(in_bam, ".bai")) && !file.exists(paste0(in_bam, ".csi"))) {
                stop("Please make sure the BAM file is indexed")
            }
            using_bam <- TRUE
            config$pipeline_parameters$do_genome_alignment <- FALSE
        }

        # setup of internal arguments which hold output files and intermediate files
        if (do_isoform_id && isoform_id_bambu) {
            isoform_gff3 <- paste(outdir, "isoform_annotated.gtf", sep = "/") # Bambu outputs GTF
        } else {
            isoform_gff3 <- paste(outdir, "isoform_annotated.gff3", sep = "/")
        }
        isoform_gff3_f <- paste(outdir, "isoform_annotated.filtered.gff3",
            sep =
                "/"
        )
        FSM_anno_out <- paste(outdir, "isoform_FSM_annotation.csv", sep = "/")
        raw_splice_isoform <- paste(outdir, "splice_raw.gff3", sep = "/")
        tss_tes_stat <- paste(outdir, "tss_tes.bedgraph", sep = "/")
        transcript_fa <- paste(outdir, "transcript_assembly.fa", sep = "/")
        transcript_fa_idx <- paste(outdir, "transcript_assembly.fa.fai",
            sep =
                "/"
        )
        tmp_bam <- paste(outdir, "tmp_align.bam", sep = "/")
        tmp_bed <- paste(outdir, "tmp_splice_anno.bed12", sep = "/")
        tmp_sam <- paste(outdir, "tmp_align.sam", sep = "/")
        genome_bam <- paste(outdir, "align2genome.bam", sep = "/")
        realign_bam <- paste(outdir, "realign2transcript.bam", sep = "/")
        tr_cnt_csv <- paste(outdir, "transcript_count.csv.gz", sep = "/")
        tr_badcov_cnt_csv <- paste(outdir, "transcript_count.bad_coverage.csv.gz",
            sep =
                "/"
        )

        cat("#### Input parameters:\n")
        print_config(config)
        cat("\tgene annotation:", annot, "\n")
        cat("\tgenome fasta:", genome_fa, "\n")
        if (using_bam) {
            cat("\tinput bam:", in_bam, "\n")
            genome_bam <- in_bam
        }
        cat("\tinput fastq:", fastq, "\n")
        cat("\toutput directory:", outdir, "\n")
        cat("\tdirectory containing minimap2:", minimap2_dir, "\n")

        # align reads to genome
        # if (!using_bam && config$pipeline_parameters$do_genome_alignment) {
        if (config$pipeline_parameters$do_genome_alignment) {
            cat("#### Aligning reads to genome using minimap2\n")

            if (config$alignment_parameters$use_junctions) {
                gff3_to_bed12(minimap2_dir, annot, tmp_bed)
            }
            minimap2_align(
                minimap2_dir,
                genome_fa,
                fastq,
                tmp_sam,
                no_flank = config$alignment_parameters$no_flank,
                bed12_junc = if (config$alignment_parameters$use_junctions) {
                    tmp_bed
                } else {
                    NULL
                },
                seed
            )
            samtools_as_bam(tmp_sam, tmp_bam)
            samtools_sort_index(tmp_bam, genome_bam)
            file.remove(tmp_sam)
            file.remove(tmp_bam)
            if (config$alignment_parameters$use_junctions) {
                file.remove(tmp_bed)
            }
        } else {
            cat("#### Skip aligning reads to genome\n")
        }

        # find isofroms
        if (isoform_id_bambu) {
            bambuAnnotations <- bambu::prepareAnnotations(annot)
            # Tmp fix: remove withr if bambu imports seqlengths properly
            # https://github.com/GoekeLab/bambu/issues/255
            bambu_out <- withr::with_package("GenomeInfoDb", bambu::bambu(reads = genome_bam, annotations = bambuAnnotations, genome = genome_fa, quant = FALSE))
            bambu::writeToGTF(bambu_out, isoform_gff3) # Does bambu_out include both novel and known isoforms ???
            # Create transcriptome assembly .fa
            gffread_cpp(genome_fa, transcript_fa, isoform_gff3) # Use XStringSet + extractTranscriptSeqs instead?
            Rsamtools::indexFa(transcript_fa)
            # Todo: convert bambu_out (GRangesList) to transcript_dict directly
            isoform_objects <- list(transcript_dict = NULL, transcript_dict_i = parse_gff_tree(isoform_gff3)$transcript_dict)
        } else {
            isoform_objects <-
                find_isoform(
                    annot,
                    genome_bam,
                    isoform_gff3,
                    tss_tes_stat,
                    genome_fa,
                    transcript_fa,
                    downsample_ratio,
                    config,
                    raw_splice_isoform
                )
        }

        # realign to transcript
        # if (!using_bam && do_read_realign) {
        if (do_read_realign) {
            cat("#### Realign to transcript using minimap2\n")
            minimap2_tr_align(minimap2_dir, transcript_fa, fastq, tmp_sam, seed)
            samtools_as_bam(tmp_sam, tmp_bam)
            samtools_sort_index(tmp_bam, realign_bam)
            file.remove(tmp_sam)
            file.remove(tmp_bam)
        } else {
            cat("#### Skip read realignment\n")
        }

        # quantification
        if (config$pipeline_parameters$do_transcript_quantification) {
            cat("#### Generating transcript count matrix\n")
            parse_realign <-
                parse_realigned_bam(
                    realign_bam,
                    transcript_fa_idx,
                    config$isoform_parameters$Min_sup_cnt,
                    config$transcript_counting$min_tr_coverage,
                    config$transcript_counting$min_read_coverage
                )
            tr_cnt <- wrt_tr_to_csv(
                parse_realign$bc_tr_count_dict,
                isoform_objects$transcript_dict_i,
                tr_cnt_csv,
                isoform_objects$transcript_dict,
                config$global_parameters$has_UMI
            )
            wrt_tr_to_csv(
                parse_realign$bc_tr_badcov_count_dict,
                isoform_objects$transcript_dict_i,
                tr_badcov_cnt_csv,
                isoform_objects$transcript_dict,
                config$global_parameters$has_UMI
            )
            annotate_filter_gff(
                isoform_gff3,
                annot,
                isoform_gff3_f,
                FSM_anno_out,
                tr_cnt,
                config$isoform_parameters$Min_sup_cnt
            )
        } else {
            cat("#### Skip transcript quantification\n")
        }

        return(
            list(
                "annot" = annot,
                "genome_fa" = genome_fa,
                "counts" = tr_cnt_csv,
                "isoform_annotated" = isoform_gff3_f,
                "transcript_assembly" = transcript_fa,
                "config" = config_file,
                "align_bam" = genome_bam,
                "realign2transcript" = realign_bam,
                "tss_tes" = tss_tes_stat,
                "outdir" = outdir
            )
        )
    }

check_arguments <-
    function(annot,
             fastq,
             in_bam,
             outdir,
             genome_fa,
             minimap2_dir,
             downsample_ratio,
             config_file,
             do_genome_align,
             do_isoform_id = TRUE,
             isoform_id_bambu = FALSE,
             do_read_realign,
             do_transcript_quanti,
             gen_raw_isoform,
             has_UMI,
             MAX_DIST,
             MAX_TS_DIST,
             MAX_SPLICE_MATCH_DIST,
             min_fl_exon_len,
             Max_site_per_splice,
             Min_sup_cnt,
             Min_cnt_pct,
             Min_sup_pct,
             strand_specific,
             remove_incomp_reads,
             use_junctions,
             no_flank,
             use_annotation,
             min_tr_coverage,
             min_read_coverage) {
        if (!dir.exists(outdir)) {
            cat("Output directory does not exists: one is being created\n")
            dir.create(outdir)
            print(outdir)
        }

        if (!do_isoform_id) {
            stop(
                "Isoform Identification is required for FLAMES execution. Change this value in the configuration file or \
        set the argument as TRUE. If isoform identification is not required, you can manually execute the pipeline by following \
        the vignette"
            )
        }

        if (is.null(config_file)) {
            config_file <- create_config(
                outdir,
                do_genome_align,
                do_isoform_id,
                do_read_realign,
                do_transcript_quanti,
                gen_raw_isoform,
                has_UMI,
                MAX_DIST,
                MAX_TS_DIST,
                MAX_SPLICE_MATCH_DIST,
                min_fl_exon_len,
                Max_site_per_splice,
                Min_sup_cnt,
                Min_cnt_pct,
                Min_sup_pct,
                strand_specific,
                remove_incomp_reads,
                use_junctions,
                no_flank,
                use_annotation,
                min_tr_coverage,
                min_read_coverage
            )
        }

        # argument verificiation
        if (downsample_ratio > 1 || downsample_ratio <= 0) {
            stop("downsample_ratio should be between 0 and 1")
        }
        if (!is.null(fastq) &&
            any(!file.exists(fastq))) {
            stop(paste0("Make sure ", fastq, " exists."))
        }
        if (!file.exists(annot)) {
            stop(paste0("Make sure ", annot, " exists."))
        }
        if (!file.exists(genome_fa)) {
            stop(paste0("Make sure ", genome_fa, " exists."))
        }

        if (!is.null(in_bam)) {
            if (any(!file.exists(in_bam))) {
                stop("Make sure in_bam exists")
            }
        }

        if (do_genome_align || do_read_realign) {
            if (!minimap2_check_callable(minimap2_dir)) {
                stop(paste0("minimap2 is not available from the given directory: ", minimap2_dir))
            }
        }

        if (isoform_id_bambu) {
            if (tail(stringr::str_split(annot, "\\.")[[1]], n = 1) != "gtf") {
                stop("Bambu requires GTF format for annotation file.\n")
            }
        }

        return(list(config = config_file))
    }