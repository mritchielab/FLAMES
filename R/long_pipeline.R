#' Generic FLAMES pipeline
#'
#' Generic implementation of the flames pipeline. Used for both bulk reads and
#' single cell reads.
#'
#' @inheritParams sc_long_pipeline
#' @param genome_bam optional BAM file which replaces fastq directory argument. This skips the genome alignment and
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
    function(annotation,
             fastq,
             genome_bam,
             outdir,
             genome_fa,
             minimap2_dir,
             config) {
        cat("Running FLAMES pipeline...\n")

        using_bam <- FALSE
        if (!is.null(genome_bam)) {
            if (!file.exists(paste0(genome_bam, ".bai")) && !file.exists(paste0(genome_bam, ".csi"))) {
                stop("Please make sure the BAM file is indexed")
            }
            using_bam <- TRUE
            config$pipeline_parameters$do_genome_alignment <- FALSE
        }

        cat("#### Input parameters:\n")
        cat(jsonlite::toJSON(config, pretty = TRUE), "\n")
        cat("gene annotation:", annotation, "\n")
        cat("genome fasta:", genome_fa, "\n")
        if (using_bam) {
            cat("input bam:", genome_bam, "\n")
        } else {
            genome_bam <- file.path(outdir, "align2genome.bam")
        }
        cat("input fastq:", fastq, "\n")
        cat("output directory:", outdir, "\n")
        cat("directory containing minimap2:", minimap2_dir, "\n")

        # align reads to genome
        # if (!using_bam && config$pipeline_parameters$do_genome_alignment) {
        if (config$pipeline_parameters$do_genome_alignment) {
            cat("#### Aligning reads to genome using minimap2\n")
            # minimap2_align <- function(config, fa_file, fq_in, annotation, outdir, minimap2_dir, threads = NULL)
            minimap2_align(
                config,
                genome_fa,
                fastq,
                annotation,
                outdir,
                minimap2_dir,
                prefix = NULL,
                threads = 12
            )
        } else {
            cat("#### Skip aligning reads to genome\n")
        }

        # find isofroms
        if (config$pipeline_parameters$do_isoform_identification) {
            find_isoform(annotation, genome_fa, genome_bam, outdir, config)
        }

        # realign to transcript
        if (config$pipeline_parameters$do_read_realignment) {
            cat("#### Realign to transcript using minimap2\n")
            minimap2_realign(config, fastq, outdir, minimap2_dir, prefix = NULL, threads = 12)
        } else {
            cat("#### Skip read realignment\n")
        }

        # quantification
        if (config$pipeline_parameters$do_transcript_quantification) {
            cat("#### Generating transcript count matrix\n")
            flames_quantify(annotation = annotation, outdir = outdir, config = config)
        } else {
            cat("#### Skip transcript quantification\n")
        }

        return(
            list(
                "annotation" = annotation,
                "genome_fa" = genome_fa,
                "counts" = file.path(outdir, "transcript_count.csv.gz"),
                "isoform_annotated" = file.path(outdir, "isoform_annotated.filtered.gff3"),
                "transcript_assembly" = file.path(outdir, "transcript_assembly.fa"),
                "align_bam" = genome_bam,
                "realign2transcript" = file.path(outdir, "realign2transcript.bam"),
                "tss_tes" = file.path(outdir, "tss_tes.bedgraph"),
                "outdir" = outdir
            )
        )
    }

#' @importFrom Matrix tail
#' @importFrom stringr str_split
#' @importFrom jsonlite fromJSON
check_arguments <-
    function(annotation,
             fastq,
             genome_bam,
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
        if (!file.exists(annotation)) {
            stop(paste0("Make sure ", annotation, " exists."))
        }
        if (!file.exists(genome_fa)) {
            stop(paste0("Make sure ", genome_fa, " exists."))
        }

        if (!is.null(genome_bam)) {
            if (any(!file.exists(genome_bam))) {
                stop("Make sure genome_bam exists")
            }
        }

        if (config$pipeline_parameters$do_genome_alignment || config$pipeline_parameters$do_read_realignment) {
            minimap2_dir <- locate_minimap2_dir(minimap2_dir = minimap2_dir)
        }

        if (config$pipeline_parameters$bambu_isoform_identification) {
            if (Matrix::tail(stringr::str_split(annotation, "\\.")[[1]], n = 1) != "gtf") {
                stop("Bambu requires GTF format for annotation file.\n")
            }
        }

        return(list(config = config, minimap2_dir = minimap2_dir))
    }
