#' Pipeline for multi-sample long-read scRNA-seq data
#'
#' @md
#'
#' @description
#' Semi-supervised isofrom detection and annotation for long read data. This variant is
#' meant for multi-sample scRNA-seq data. Specific parameters can be configured in
#' the config file (see \code{\link{create_config}}), input files are specified via
#' arguments.
#'
#' @inherit SingleCellPipeline details
#' @inheritParams SingleCellPipeline
#' @param fastq A named vector of fastq file (or folder) paths. Each element of the vector
#'   will be treated as a sample. The names of the vector will be used as the sample names.
#'   If not named, the sample names will be generated from the file names.
#'
#' @return A \code{FLAMES.MultiSampleSCPipeline} object. The pipeline can be run using
#'   the \code{\link{run_FLAMES}} function. The resulting list of SingleCellExperiment
#'   objects can be accessed using the \code{experiment} method.
#'
#' @seealso
#' \code{\link{SingleCellPipeline}} for single-sample long data and more details on the
#' pipeline output,
#' \code{\link{create_config}} for creating a configuration file,
#' \code{\link{BulkPipeline}} for bulk long data.
#'
#' @examples
#' reads <- ShortRead::readFastq(
#'   system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
#' )
#' outdir <- tempfile()
#' dir.create(outdir)
#' dir.create(file.path(outdir, "fastq"))
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' genome_fa <- file.path(outdir, "rps24.fa")
#' R.utils::gunzip(
#'   filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
#'   destname = bc_allow, remove = FALSE
#' )
#' R.utils::gunzip(
#'   filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
#'   destname = genome_fa, remove = FALSE
#' )
#' ShortRead::writeFastq(reads[1:100],
#'   file.path(outdir, "fastq/sample1.fq.gz"), mode = "w", full = FALSE)
#' reads <- reads[-(1:100)]
#' ShortRead::writeFastq(reads[1:100],
#'   file.path(outdir, "fastq/sample2.fq.gz"), mode = "w", full = FALSE)
#' reads <- reads[-(1:100)]
#' ShortRead::writeFastq(reads,
#'   file.path(outdir, "fastq/sample3.fq.gz"), mode = "w", full = FALSE)
#' ppl <- MultiSampleSCPipeline(
#'   config_file = create_config(
#'     outdir,
#'     pipeline_parameters.demultiplexer = "flexiplex",
#'     pipeline_parameters.threads = 1,
#'     alignment_parameters.no_flank = TRUE
#'   ),
#'   outdir = outdir,
#'   fastq = c("sampleA" = file.path(outdir, "fastq"),
#'     "sample1" = file.path(outdir, "fastq", "sample1.fq.gz"),
#'     "sample2" = file.path(outdir, "fastq", "sample2.fq.gz"),
#'     "sample3" = file.path(outdir, "fastq", "sample3.fq.gz")),
#'   annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
#'   genome_fa = genome_fa,
#'   barcodes_file = rep(bc_allow, 4)
#' )
#' ppl <- run_FLAMES(ppl)
#' experiment(ppl)
#' @export
MultiSampleSCPipeline <- function(
  config_file, outdir, fastq, annotation, genome_fa, genome_mmi,
  minimap2, samtools, barcodes_file, expect_cell_number, controllers
) {
  pipeline <- new("FLAMES.MultiSampleSCPipeline")
  config <- check_arguments(annotation, fastq, genome_bam = NULL, outdir, genome_fa, config_file)$config

  if (!dir.exists(outdir)) {
    dir.create(outdir)
    message(sprintf("Output directory (%s) did not exist, created.", outdir))
  }

  if (length(fastq) == 1) {
    if (file_test("-f", fastq)) {
      stop("Only one fastq file provided, did you meant to used the single-sample pipeline (FLAMES::sc_long_pipeline) ?")
    }
    fastq <- list.files(fastq, pattern = "\\.(fastq|fq)(\\.gz)?$", full.names = TRUE)
    if (length(fastq) <= 1) {
      stop(length(fastq), " .fq or .fastq file(s) found\n")
    }
  } else if (!all(file.exists(fastq))) {
    stop("fastq must be a valid path to a folder or a FASTQ file")
  }
  if (is.null(names(fastq))) {
    names(fastq) <- gsub("\\.(fastq|fq)(\\.gz)?$", "", basename(fastq))
  }
  names(fastq) <- make.unique(names(fastq), sep = "_")

  steps <- c(
    "barcode_demultiplex", "genome_alignment", "gene_quantification",
    "isoform_identification", "read_realignment", "transcript_quantification"
  )
  steps <- config$pipeline_parameters[paste0("do_", steps)] |>
    unlist() |>
    setNames(steps)
  message("Configured steps: \n", paste0("\t", names(steps), ": ", steps, collapse = "\n"))

  # assign slots
  ## inputs
  pipeline@config <- config
  pipeline@outdir <- outdir
  pipeline@fastq <- fastq
  pipeline@annotation <- annotation
  pipeline@genome_fa <- genome_fa
  if (!missing(genome_mmi) && is.character(genome_mmi)) {
    pipeline@genome_mmi <- genome_mmi
  }
  if (!missing(barcodes_file) && is.character(barcodes_file)) {
    if (length(barcodes_file) == 1) {
      barcodes_file <- rep(barcodes_file, length(fastq))
    } else if (length(barcodes_file) != length(fastq)) {
      stop(
        sprintf(
          "Number of barcodes_file (%d) does not match number of fastq files (%d)",
          length(barcodes_file), length(fastq)
        )
      )
    }
    pipeline@barcodes_file <- barcodes_file
  }
  if (!missing(expect_cell_number) && is.numeric(expect_cell_number)) {
    if (length(expect_cell_number) == 1) {
      expect_cell_number <- rep(expect_cell_number, length(fastq))
    } else if (length(expect_cell_number) != length(fastq)) {
      stop(
        sprintf(
          "Number of expect_cell_number (%d) does not match number of fastq files (%d)",
          length(expect_cell_number), length(fastq)
        )
      )
    }
    pipeline@expect_cell_number <- expect_cell_number
  }

  if (steps["barcode_demultiplex"] &&
    config$pipeline_parameters$demultiplexer == "BLAZE" &&
    any(is.na(pipeline@expect_cell_number))) {
    stop("demultiplexer set to 'BLAZE', which requires 'expect_cell_number' to be provided.")
  }

  ## outputs
  # metadata
  pipeline@bed <- file.path(outdir, "reference.bed")
  pipeline@genome_bam <- file.path(outdir, paste0(names(fastq), "_", "align2genome.bam"))
  pipeline@transcriptome_bam <- file.path(outdir, paste0(names(fastq), "_", "realign2transcript.bam"))
  pipeline@transcriptome_assembly <- file.path(outdir, "transcript_assembly.fa")
  pipeline@demultiplexed_fastq <- file.path(outdir, paste0(names(fastq), "_matched_reads.fastq.gz"))
  pipeline@deduped_fastq <- file.path(outdir, paste0(names(fastq), "_matched_reads_dedup.fastq.gz"))
  pipeline@experiment <- file.path(outdir, "experiment.rds")

  ## binaries
  if (missing(minimap2) || !is.character(minimap2)) {
    minimap2 <- find_bin("minimap2")
    if (is.na(minimap2)) {
      stop("minimap2 not found, please make sure it is installed and provide its path as the minimap2 argument")
    }
  }
  pipeline@minimap2 <- minimap2
  if (missing(samtools) || !is.character(samtools)) {
    samtools <- find_bin("samtools")
  }
  if (is.na(samtools)) {
    message("samtools not found, will use Rsamtools package instead")
  }
  pipeline@samtools <- samtools
  ##

  ## pipeline state
  pipeline@steps <- steps
  pipeline@completed_steps <- setNames(
    rep(FALSE, length(steps)), names(steps)
  )

  if (!missing(controllers)) {
    pipeline@controllers <- normalize_controllers(controllers, names(steps))
  }

  return(pipeline)
}

setMethod("barcode_demultiplex", "FLAMES.MultiSampleSCPipeline", function(pipeline) {
  using_controllers <- FALSE
  if ("barcode_demultiplex" %in% names(pipeline@controllers)) {
    using_controllers <- TRUE
    controller <- pipeline@controllers[["barcode_demultiplex"]]
  } else if ("default" %in% names(pipeline@controllers)) {
    using_controllers <- TRUE
    controller <- pipeline@controllers[["default"]]
  }

  demux_tool <- pipeline@config$pipeline_parameters$demultiplexer

  if (demux_tool == "BLAZE") {
    message("Using BLAZE for barcode demultiplexing.")
    if (any(pipeline@barcodes_file != "")) {
        warning("`barcodes_file` is provided for one or more samples but demultiplexer is set to 'BLAZE'. The file(s) will be ignored.")
    }
    if (any(is.na(pipeline@expect_cell_number))) {
        stop("demultiplexer is 'BLAZE', so 'expect_cell_number' must be provided for all samples.")
    }

    if (using_controllers) {
      if (controller$started()) {
        tryCatch({
          controller$terminate()
        }, error = function(e) {
          warning(sprintf("Error terminating controller: %s", e$message))
        })
      }
      controller$start()
      controller$map(
        command = blaze(
          expect_cell_number[i], fastq[i],
          additional_args = config$additional_arguments$BLAZE,
          "output-fastq" = demultiplexed_fastq[i],
          "threads" = config$pipeline_parameters$threads,
          "max-edit-distance" = config$barcode_parameters$max_bc_editdistance,
          "overwrite" = TRUE
        ),
        iterate = list(i = seq_along(pipeline@fastq)),
        data = list(
          expect_cell_number = pipeline@expect_cell_number,
          fastq = pipeline@fastq,
          outdir = pipeline@outdir,
          demultiplexed_fastq = pipeline@demultiplexed_fastq,
          config = pipeline@config,
          blaze = FLAMES:::blaze
        ),
        error = "stop"
      )
      controller$terminate()
    } else {
      for (i in seq_along(pipeline@fastq)) {
        message(sprintf("Demultiplexing sample %s...", names(pipeline@fastq)[i]))
        blaze(
          pipeline@expect_cell_number[i], pipeline@fastq[i],
          additional_args = pipeline@config$additional_arguments$BLAZE,
          "output-fastq" = pipeline@demultiplexed_fastq[i],
          "threads" = pipeline@config$pipeline_parameters$threads,
          "max-edit-distance" = pipeline@config$barcode_parameters$max_bc_editdistance,
          "overwrite" = TRUE
        )
      }
    }
  } else if (demux_tool == "flexiplex") {
    message("Using flexiplex for barcode demultiplexing.")
    if (any(!is.na(pipeline@expect_cell_number))) {
        warning("'expect_cell_number' is provided for one or more samples but is not used by flexiplex.")
    }

    if (using_controllers) {
      if (controller$started()) {
        tryCatch({
          controller$terminate()
        }, error = function(e) {
          warning(sprintf("Error terminating controller: %s", e$message))
        })
      }
      controller$start()
      crew_result <- controller$map(
        command = find_barcode(
          fastq = fastq[i],
          barcodes_file = barcodes_file[i],
          stats_out = file.path(
            outdir,
            paste0(names(fastq), "_matched_barcode_stat.tsv.gz")
          )[i],
          reads_out = demultiplexed_fastq[i],
          pattern = setNames(
            as.character(config$barcode_parameters$pattern),
            names(config$barcode_parameters$pattern)
          ),
          TSO_seq = config$barcode_parameters$TSO_seq,
          TSO_prime = config$barcode_parameters$TSO_prime,
          cutadapt_minimum_length = config$barcode_parameters$cutadapt_minimum_length,
          full_length_only = config$barcode_parameters$full_length_only,
          max_bc_editdistance = config$barcode_parameters$max_bc_editdistance,
          max_flank_editdistance = config$barcode_parameters$max_flank_editdistance,
          strand = config$barcode_parameters$strand,
          threads = config$pipeline_parameters$threads
        ),
        iterate = list(i = seq_along(pipeline@fastq)),
        data = list(
          fastq = pipeline@fastq,
          barcodes_file = pipeline@barcodes_file,
          outdir = pipeline@outdir,
          demultiplexed_fastq = pipeline@demultiplexed_fastq,
          config = pipeline@config,
          find_barcode = FLAMES:::find_barcode
        ),
        error = "stop"
      )
      controller$terminate()
      res <- crew_result$result
    } else {
      res <- lapply(seq_along(pipeline@fastq), function(i) {
        find_barcode(
          fastq = pipeline@fastq[i],
          barcodes_file = pipeline@barcodes_file[i],
          stats_out = file.path(
            pipeline@outdir,
            paste0(names(pipeline@fastq)[i], "_matched_barcode_stat.tsv.gz")
          ),
          reads_out = pipeline@demultiplexed_fastq[i],
          pattern = setNames(
            as.character(pipeline@config$barcode_parameters$pattern),
            names(pipeline@config$barcode_parameters$pattern)
          ),
          TSO_seq = pipeline@config$barcode_parameters$TSO_seq,
          TSO_prime = pipeline@config$barcode_parameters$TSO_prime,
          cutadapt_minimum_length = pipeline@config$barcode_parameters$cutadapt_minimum_length,
          full_length_only = pipeline@config$barcode_parameters$full_length_only,
          max_bc_editdistance = pipeline@config$barcode_parameters$max_bc_editdistance,
          max_flank_editdistance = pipeline@config$barcode_parameters$max_flank_editdistance,
          strand = pipeline@config$barcode_parameters$strand,
          threads = pipeline@config$pipeline_parameters$threads
        )
      })
      names(res) <- names(pipeline@fastq)
    }
    pipeline@metadata$barcode_demultiplex <- res
  } else {
    stop("Unknown demultiplexer specified in config: ", demux_tool)
  }
  return(pipeline)
})

#' Pipeline for Multi-sample Single Cell Data (deprecated)
#'
#' @description This function is deprecated. Please use \code{\link{MultiSampleSCPipeline}}.
#'
#' @param annotation The file path to the annotation file in GFF3 format
#' @param fastqs The file path to input fastq file
#' @param outdir The path to directory to store all output files.
#' @param genome_fa The file path to genome fasta file.
#' @param minimap2 Path to minimap2, optional.
#' @param barcodes_file The file with expected cell barcodes, with each barcode on a new line.
#' @param expect_cell_numbers The expected number of cells in the sample. This is used if
#'   \code{barcodes_file} is not provided. See \code{BLAZE} for more details.
#' @param config_file File path to the JSON configuration file.
#'
#' @return A list of \code{SingleCellExperiment} objects, one for each sample.
#'
#' @seealso
#' \code{\link{MultiSampleSCPipeline}} for the new pipeline interface,
#' \code{\link{SingleCellPipeline}} for single-sample pipeline,
#' \code{\link{BulkPipeline}} for bulk long data.
#'
#' @importFrom utils file_test
#'
#' @examples
#' reads <- ShortRead::readFastq(
#'   system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
#' )
#' outdir <- tempfile()
#' dir.create(outdir)
#' dir.create(file.path(outdir, "fastq"))
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' genome_fa <- file.path(outdir, "rps24.fa")
#' R.utils::gunzip(
#'   filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
#'   destname = bc_allow, remove = FALSE
#' )
#' R.utils::gunzip(
#'   filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
#'   destname = genome_fa, remove = FALSE
#' )
#' ShortRead::writeFastq(reads[1:100],
#'   file.path(outdir, "fastq/sample1.fq.gz"), mode = "w", full = FALSE)
#' reads <- reads[-(1:100)]
#' ShortRead::writeFastq(reads[1:100],
#'   file.path(outdir, "fastq/sample2.fq.gz"), mode = "w", full = FALSE)
#' reads <- reads[-(1:100)]
#' ShortRead::writeFastq(reads,
#'   file.path(outdir, "fastq/sample3.fq.gz"), mode = "w", full = FALSE)
#'
#' sce_list <- FLAMES::sc_long_multisample_pipeline(
#'   annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
#'   fastqs = c("sampleA" = file.path(outdir, "fastq"),
#'     "sample1" = file.path(outdir, "fastq", "sample1.fq.gz"),
#'     "sample2" = file.path(outdir, "fastq", "sample2.fq.gz"),
#'     "sample3" = file.path(outdir, "fastq", "sample3.fq.gz")),
#'   outdir = outdir,
#'   genome_fa = genome_fa,
#'   barcodes_file = rep(bc_allow, 4),
#'   config_file = create_config(
#'     outdir,
#'     pipeline_parameters.demultiplexer = "flexiplex"
#'   )
#' )
#'
#' @export
sc_long_multisample_pipeline <- function(annotation, fastqs, outdir, genome_fa,
    minimap2 = NULL, barcodes_file = NULL,
    expect_cell_numbers = NULL, config_file = NULL) {
  message("sc_long_multisample_pipeline is deprecated, please use MultiSampleSCPipeline instead.")
  pipeline <- MultiSampleSCPipeline(
    config_file = config_file,
    outdir = outdir,
    fastq = fastqs,
    annotation = annotation,
    genome_fa = genome_fa,
    minimap2 = minimap2,
    barcodes_file = barcodes_file,
    expect_cell_number = expect_cell_numbers
  )
  pipeline <- run_FLAMES(pipeline)
  saveRDS(pipeline, file.path(outdir, "pipeline.rds"))
  message("Pipeline saved to ", file.path(outdir, "pipeline.rds"))
  if (length(pipeline@last_error) == 0) {
    return(experiment(pipeline))
  } else {
    warning("Returning pipeline object instead of experiment due to errors.")
    message("You can resume the pipeline after resolving the errors with resume_FLAMES(pipeline)")
    return(pipeline)
  }
}
