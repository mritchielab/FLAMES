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
#' @importFrom jsonlite fromJSON toJSON
generic_long_pipeline <-
    function(annot,
             fastq,
             in_bam,
             outdir,
             genome_fa,
             minimap2_dir,
             config) {
        cat("Running FLAMES pipeline...\n")

        using_bam <- FALSE
        if (!is.null(in_bam)) {
            if (!file.exists(paste0(in_bam, ".bai")) && !file.exists(paste0(in_bam, ".csi"))) {
                stop("Please make sure the BAM file is indexed")
            }
            using_bam <- TRUE
            config$pipeline_parameters$do_genome_alignment <- FALSE
        }

        # setup of internal arguments which hold output files and intermediate files
        if (config$pipeline_parameters$do_isoform_identification && config$pipeline_parameters$bambu_isoform_identification) {
            isoform_gff3 <- paste(outdir, "isoform_annotated.gtf", sep = "/") # Bambu outputs GTF
        } else {
            isoform_gff3 <- paste(outdir, "isoform_annotated.gff3", sep = "/")
        }
        isoform_gff3_f <- paste(outdir, "isoform_annotated.filtered.gff3",
            sep =
                "/"
        )
        transcript_fa <- paste(outdir, "transcript_assembly.fa", sep = "/")
        transcript_fa_idx <- paste(outdir, "transcript_assembly.fa.fai",
            sep =
                "/"
        )
        genome_bam <- paste(outdir, "align2genome.bam", sep = "/")
        realign_bam <- paste(outdir, "realign2transcript.bam", sep = "/")
        tr_cnt_csv <- paste(outdir, "transcript_count.csv.gz", sep = "/")
        tr_badcov_cnt_csv <- paste(outdir, "transcript_count.bad_coverage.csv.gz",
            sep =
                "/"
        )

        cat("#### Input parameters:\n")
        cat(jsonlite::toJSON(config, pretty = TRUE), "\n")
        cat("gene annotation:", annot, "\n")
        cat("genome fasta:", genome_fa, "\n")
        if (using_bam) {
            cat("input bam:", in_bam, "\n")
            genome_bam <- in_bam
        }
        cat("input fastq:", fastq, "\n")
        cat("output directory:", outdir, "\n")
        cat("directory containing minimap2:", minimap2_dir, "\n")

        # align reads to genome
        # if (!using_bam && config$pipeline_parameters$do_genome_alignment) {
        if (config$pipeline_parameters$do_genome_alignment) {
            cat("#### Aligning reads to genome using minimap2\n")
            # minimap2_align <- function(config, fa_file, fq_in, annot, outdir, minimap2_dir, threads = NULL)
            minimap2_align(
                config,
                genome_fa,
                fastq,
                annot,
                outdir,
                minimap2_dir,
                prefix = NULL,
                threads = 12
            )
        } else {
            cat("#### Skip aligning reads to genome\n")
        }

        # find isofroms
        if (config$pipeline_parameters$bambu_isoform_identification) {
            bambuAnnotations <- bambu::prepareAnnotations(annot)
            # Tmp fix: remove withr if bambu imports seqlengths properly
            # https://github.com/GoekeLab/bambu/issues/255
            bambu_out <- withr::with_package("GenomeInfoDb", bambu::bambu(reads = genome_bam, annotations = bambuAnnotations, genome = genome_fa, quant = FALSE))
            bambu::writeToGTF(bambu_out, isoform_gff3) # Does bambu_out include both novel and known isoforms ???

            # Create transcriptome assembly .fa
            dna_string_set <- Biostrings::readDNAStringSet(genome_fa)
            names(dna_string_set) <- gsub(" .*$", "", names(dna_string_set))
            tr_string_set <- GenomicFeatures::extractTranscriptSeqs(dna_string_set, get_GRangesList(isoform_gff3))
            Biostrings::writeXStringSet(tr_string_set, transcript_fa)

            Rsamtools::indexFa(transcript_fa)
            # Todo: convert bambu_out (GRangesList) to transcript_dict directly
            isoform_objects <- list(transcript_dict = NULL, transcript_dict_i = parse_gff_tree(isoform_gff3)$transcript_dict)
        } else {
            isoform_objects <-
                find_isoform(
                    annot,
                    genome_bam,
                    isoform_gff3,
                    file.path(outdir, "tss_tes.bedgraph"),
                    genome_fa,
                    transcript_fa,
                    config$isoform_parameters$downsample_ratio,
                    config,
                    file.path(outdir, "splice_raw.gff3")
                )
        }

        # realign to transcript
        if (config$pipeline_parameters$do_read_realignment) {
            cat("#### Realign to transcript using minimap2\n")
            minimap2_realign(config, transcript_fa, fastq, outdir, minimap2_dir, prefix = NULL, threads = 12)
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
                file.path(outdir, "isoform_FSM_annotation.csv"),
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
                "align_bam" = genome_bam,
                "realign2transcript" = realign_bam,
                "tss_tes" = file.path(outdir, "tss_tes.bedgraph"),
                "outdir" = outdir
            )
        )
    }

#' @importFrom Matrix tail
#' @importFrom stringr str_split
#' @importFrom jsonlite fromJSON
check_arguments <-
    function(annot,
             fastq,
             in_bam,
             outdir,
             genome_fa,
             minimap2_dir,
             config_file) {
        if (!dir.exists(outdir)) {
            cat("Output directory does not exists: one is being created\n")
            dir.create(outdir)
            print(outdir)
        }

        if (is.null(config_file)) {
            cat("No config file provided, creating a default config in", outdir, "\n")
            config_file <- create_config(outdir)
        }

        # argument verificiation
        config <- jsonlite::fromJSON(config_file)

        if (config$isoform_parameters$downsample_ratio > 1 || config$isoform_parameters$downsample_ratio <= 0) {
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

        if (config$pipeline_parameters$do_genome_alignment || config$pipeline_parameters$do_read_realignment) {
            minimap2_dir <- locate_minimap2_dir(minimap2_dir = minimap2_dir)
        }

        if (config$pipeline_parameters$bambu_isoform_identification) {
            if (Matrix::tail(stringr::str_split(annot, "\\.")[[1]], n = 1) != "gtf") {
                stop("Bambu requires GTF format for annotation file.\n")
            }
        }

        return(list(config = config, minimap2_dir = minimap2_dir))
    }

get_GRangesList <- function(file, gene = NULL) {
    isoform_gff <- rtracklayer::import(file, feature.type = c("exon", "utr"))
    if (is.null("gene")) {
        isoform_gff <- isoform_gff[isoform_gff$gene_id == gene]
    }
    isoform_gff <- S4Vectors::split(isoform_gff, isoform_gff$transcript_id)
    isoform_gff
}
