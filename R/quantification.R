#' @importFrom reticulate import_from_path dict
#' @importFrom basilisk basiliskRun
parse_realigned_bam <- function(bam_in, fa_idx_f, min_sup_reads,
    min_tr_coverage, min_read_coverage, bc_file) {
  ret <- basiliskRun(env = flames_env,
    fun = function(bam_in, fa_idx_f, min_sup_reads,
        min_tr_coverage, min_read_coverage, bc_file) {
      python_path <- system.file("python", package = "FLAMES")
      count <- reticulate::import_from_path("count_tr", python_path)
      ret <- count$parse_realigned_bam(
        bam_in, fa_idx_f, min_sup_reads,
        min_tr_coverage, min_read_coverage, bc_file
      )
      names(ret) <- c("bc_tr_count_dict", "bc_tr_badcov_count_dict", "tr_kept")
      ret
    },
    bam_in = bam_in, fa_idx_f = fa_idx_f, min_sup_reads = min_sup_reads,
    min_tr_coverage = min_tr_coverage, min_read_coverage = min_read_coverage,
    bc_file = bc_file
  )
  ret
}

#' @importFrom reticulate import_from_path
#' @importFrom basilisk basiliskRun
#' @importFrom future plan
wrt_tr_to_csv <- function(bc_tr_count_dict, transcript_dict, csv_f,
    transcript_dict_ref = NULL, has_UMI = TRUE) {
  future::plan(future::multisession)
  basiliskRun(env = flames_env,
    fun = function(bc_tr_count_dict, transcript_dict,
        csv_f, transcript_dict_ref, has_UMI) {
      python_path <- system.file("python", package = "FLAMES")
      count <- reticulate::import_from_path("count_tr", python_path)
      count$wrt_tr_to_csv(bc_tr_count_dict, transcript_dict, csv_f, transcript_dict_ref, has_UMI)
    },
    bc_tr_count_dict = bc_tr_count_dict, transcript_dict = transcript_dict,
    csv_f = csv_f, transcript_dict_ref = transcript_dict_ref, has_UMI = has_UMI
  )
}


#' Gene quantification
#' @description Calculate the per gene UMI count matrix by parsing the genome alignment file.
#'
#' @details
#' After the genome alignment step (\code{do_genome_align}), the alignment file will be parsed to
#' generate the per gene UMI count matrix. For each gene in the annotation file, the number of
#' reads whose mapped ranges overlap with the gene's genome coordinates will be assigned to the
#' gene. For reads can be assigned to multiple gene, the read will be assigned to the gene with
#' the highest number of overlapping nucleotides. If the read can be assigned to multiple genes
#' with the same number of overlapping nucleotides, the read will be not be assigned.
#'
#' After the read-to-gene assignment, the per gene UMI count matrix will be generated.
#' Specifically, for each gene, the reads with similar mapping coordinates of transcript
#' termination sites (TTS, i.e. the end of the the read with a polyT or polyA) will be grouped
#' together. UMIs of reads in the same group will be collapsed to generate the UMI counts for each
#' gene.
#'
#' Finally, a new fastq file with deduplicated reads by keeping the longest read in each UMI.
#'
#' @param annotation The file path to the annotation file in GFF3 format
#' @param outdir The path to directory to store all output files.
#' @param infq The input FASTQ file.
#' @param n_process The number of processes to use for parallelization.
#' @param pipeline The pipeline type as a character string, either \code{sc_single_sample} (single-cell, single-sample),
#' @param samples A vector of sample names, default to the file names of input fastq files,
#' or folder names if \code{fastqs} is a vector of folders.
#' \code{bulk} (bulk, single or multi-sample), or \code{sc_multi_sample} (single-cell, multiple samples)
#' @return The count matrix will be saved in the output folder as \code{transcript_count.csv.gz}.
#' @importFrom reticulate import_from_path dict
#' @importFrom basilisk basiliskRun
quantify_gene <- function(annotation, outdir, infq, n_process, pipeline = "sc_single_sample", samples = NULL) {
  cat(format(Sys.time(), "%X %a %b %d %Y"), "quantify genes \n")

  if (grepl("\\.gff3?(\\.gz)?$", annotation)) {
    warning("Annotation in GFF format may cause errors. Please consider using GTF formats.\n")
  }

  genome_bam <- list.files(outdir)[grepl("_?align2genome\\.bam$", list.files(outdir))]
  cat("Found genome alignment file(s): ")
  cat(paste0("\t", paste(genome_bam, collapse = "\n\t"), "\n"))

  if (length(genome_bam) != 1 && grepl("single_sample", pipeline)) {
    stop("Incorrect number of genome alignment files found.\n")
  }

  tryCatch(
    {
      basiliskRun(
        env = flames_env, fun = function(annotation, outdir, pipeline, n_process, infq, samples) {
          python_path <- system.file("python", package = "FLAMES")
          count <- reticulate::import_from_path("count_gene", python_path)
          count$quantification(annotation, outdir, pipeline, n_process, infq = infq, sample_names = samples)
        },
        annotation = annotation,
        outdir = outdir,
        pipeline = pipeline,
        n_process = n_process,
        infq = infq,
        samples = samples
      )
    },
    error = function(e) {
      # Capture the Python error using py_last_error()
      py_error <- reticulate::py_last_error()
      py_error_message <- py_error$message
      stop("Error when quantifying genes:\n", py_error_message)
    }
  )
}

#' FLAMES Transcript quantification
#' @description Calculate the transcript count matrix by parsing the re-alignment file.
#' @param annotation The file path to the annotation file in GFF3 format
#' @param outdir The path to directory to store all output files.
#' @param config Parsed FLAMES configurations.
#' @param pipeline The pipeline type as a character string, either \code{sc_single_sample} (single-cell, single-sample),
#' @param samples A vector of sample names, required for \code{sc_multi_sample} pipeline.
#' \code{bulk} (bulk, single or multi-sample), or \code{sc_multi_sample} (single-cell, multiple samples)
#' @return The count matrix will be saved in the output folder as \code{transcript_count.csv.gz}.
#' @importFrom basilisk basiliskRun
#' @importFrom reticulate import_from_path dict
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <- "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq1", paste(file_url, "fastq/sample1.fastq.gz", sep = "/")))]]
#' genome_fa <- bfc[[names(BiocFileCache::bfcadd(bfc, "genome.fa", paste(file_url, "SIRV_isoforms_multi-fasta_170612a.fasta", sep = "/")))]]
#' annotation <- bfc[[names(BiocFileCache::bfcadd(bfc, "annot.gtf", paste(file_url, "SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf", sep = "/")))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' fasta <- annotation_to_fasta(annotation, genome_fa, outdir)
#' config <- jsonlite::fromJSON(create_config(outdir, bambu_isoform_identification = TRUE, min_tr_coverage = 0.1, min_read_coverage = 0.1, min_sup_cnt = 1))
#' file.copy(annotation, file.path(outdir, "isoform_annotated.gtf"))
#' \dontrun{
#' if (!any(is.na(find_bin(c("minimap2", "k8"))))) {
#'   minimap2_realign(
#'     config = config, outdir = outdir,
#'     fq_in = fastq1
#'   )
#'   quantify_transcript_flames(annotation, outdir, config, pipeline = "bulk")
#' }
#' }
#' @export
quantify_transcript_flames <- function(annotation, outdir, config, pipeline = "sc_single_sample", samples) {
  cat(format(Sys.time(), "%X %a %b %d %Y"), "quantify transcripts \n")

  if (grepl("\\.gff3?(\\.gz)?$", annotation)) {
    warning("Annotation in GFF format may cause errors. Please consider using GTF formats.\n")
  }

  realign_bam <- list.files(outdir)[grepl("_?realign2transcript\\.bam$", list.files(outdir))]
  cat("Found realignment file(s): ")
  cat(paste0("\t", paste(realign_bam, collapse = "\n\t"), "\n"))

  if (length(realign_bam) != 1 && grepl("single_sample", pipeline)) {
    stop("Incorrect number of realignment files found.\n")
  }

  basiliskRun(
    env = flames_env, fun = function(config_dict, annotation, outdir, pipeline) {
      python_path <- system.file("python", package = "FLAMES")
      count <- reticulate::import_from_path("count_tr", python_path)
      count$quantification(config_dict, annotation, outdir, pipeline)
    },
    config_dict = reticulate::dict(config),
    annotation = annotation,
    outdir = outdir,
    pipeline = pipeline
  )

  if (pipeline == "sc_single_sample") {
    out_files <- list(
      "annotation" = annotation,
      # "genome_fa" = genome_fa,
      "counts" = file.path(outdir, "transcript_count.csv.gz"),
      "isoform_annotated" = file.path(outdir, "isoform_annotated.filtered.gff3"),
      "transcript_assembly" = file.path(outdir, "transcript_assembly.fa"),
      # "align_bam" = genome_bam,
      "realign2transcript" = file.path(outdir, "realign2transcript.bam"),
      # "tss_tes" = file.path(outdir, "tss_tes.bedgraph"),
      "fsm_annotation" = file.path(outdir, "isoform_FSM_annotation.csv"),
      "outdir" = outdir
    )
    load_genome_anno <- rtracklayer::import(annotation, feature.type = c("exon", "utr"))
    sce <- generate_sc_singlecell(out_files, load_genome_anno = load_genome_anno)
    return(sce)
  } else if (pipeline == "sc_multi_sample") {
    sce_list <- as.list(1:length(samples))
    names(sce_list) <- samples
    load_genome_anno <- rtracklayer::import(annotation, feature.type = c("exon", "utr"))
    for (i in 1:length(samples)) {
      out_files <- list(
        "annotation" = annotation,
        "counts" = file.path(outdir, paste0(samples[i], "_transcript_count.csv.gz")),
        "isoform_annotated" = file.path(outdir, ifelse(config$pipeline_parameters$bambu_isoform_identification, "isoform_annotated.gtf", "isoform_annotated.gff3")),
        "transcript_assembly" = file.path(outdir, "transcript_assembly.fa"),
        "realign2transcript" = file.path(outdir, paste0(samples[i], "_realign2transcript.bam")),
        "outdir" = outdir,
        "fsm_annotation" = file.path(outdir, "isoform_FSM_annotation.csv")
      )
      sce_list[[i]] <- generate_sc_singlecell(out_files, load_genome_anno = load_genome_anno)
    }
    return(sce_list)
  } else if (pipeline == "bulk") {
    out_files <- list(
      "annotation" = annotation,
      "counts" = file.path(outdir, "transcript_count.csv.gz"),
      "isoform_annotated" = file.path(outdir, "isoform_annotated.filtered.gff3"),
      "transcript_assembly" = file.path(outdir, "transcript_assembly.fa"),
      "realign2transcript" = file.path(outdir, list.files(outdir))[grepl("realign2transcript\\.bam$", list.files(outdir))],
      "outdir" = outdir,
      "fsm_annotation" = file.path(outdir, "isoform_FSM_annotation.csv")
    )
    load_genome_anno <- rtracklayer::import(annotation, feature.type = c("exon", "utr"))
    se <- generate_bulk_summarized(out_files, load_genome_anno = load_genome_anno)
    return(se)
  }
}

#' @importFrom basilisk obtainEnvironmentPath basiliskRun
run_oarfish <- function(realign_bam, outdir, threads = 1, sample, oarfish_bin, single_cell = TRUE) {
  if (missing(oarfish_bin)) {
    oarfish_bin <- find_bin("oarfish")
    stopifnot(!is.na(oarfish_bin))
  }

  if (missing(sample)) {
    sample <- "oarfish"
  }

  oarfish_status <- base::system2(
    command = oarfish_bin,
    args = c(
      switch(single_cell,
        "--single-cell"
      ),
      "--alignments", file.path(outdir, realign_bam),
      "-j", threads, "--output", file.path(outdir, sample)
    ),
  )
  if (oarfish_status != 0) {
    stop(paste0("error running oarfish:\n", oarfish_status))
  }

  return(file.path(outdir, sample))
}

#' @importFrom MatrixGenerics rowSums
#' @importFrom Matrix readMM
#' @importFrom SingleCellExperiment SingleCellExperiment
parse_oarfish_sc_output <- function(oarfish_out) {
  mtx <- t(Matrix::readMM(paste0(oarfish_out, ".count.mtx")))
  rownames(mtx) <- read.delim(paste0(oarfish_out, ".features.txt"), header = F)$V1
  # mtx <- mtx[MatrixGenerics::rowSums(mtx) > 0, ]
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mtx))
}

#' @importFrom SummarizedExperiment SummarizedExperiment
parse_oarfish_bulk_output <- function(oarfish_outs, sample_names) {
  mtx_list <- lapply(oarfish_outs, function(oarfish_out) {
    read.delim(paste0(oarfish_out, ".quant"), header = T, row.names = 1)[, "num_reads", drop = F]
  })
  mtx <- do.call(cbind, mtx_list)
  colnames(mtx) <- sample_names
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = mtx))
}

quantify_transcript_oarfish <- function(outdir, config, pipeline = "sc_single_sample", samples) {
  realign_bam <- list.files(outdir)[grepl("_?realign2transcript\\.bam$", list.files(outdir))]
  if (pipeline == "sc_single_sample") {
    oarfish_out <- run_oarfish(realign_bam, outdir, threads = config$pipeline_parameters$threads)
    return(parse_oarfish_sc_output(oarfish_out))
  } else if (pipeline == "sc_multi_sample") {
    sce_list <- as.list(1:length(samples))
    names(sce_list) <- samples
    for (i in 1:length(samples)) {
      realign_bam <- list.files(outdir)[grepl(paste0(samples[i], "_realign2transcript\\.bam$"), list.files(outdir))]
      oarfish_out <- run_oarfish(realign_bam, outdir, threads = config$pipeline_parameters$threads, sample = samples[i])
      sce_list[[i]] <- parse_oarfish_sc_output(oarfish_out)
    }
    return(sce_list)
  } else if (pipeline == "bulk") {
    oarfish_out <- rep(NA, length(samples))
    for (i in 1:length(samples)) {
      realign_bam <- list.files(outdir)[grepl(paste0(samples[i], "_realign2transcript\\.bam$"), list.files(outdir))]
      oarfish_out[i] <- run_oarfish(realign_bam, outdir, threads = config$pipeline_parameters$threads, sample = samples[i], single_cell = FALSE)
    }
    return(parse_oarfish_bulk_output(oarfish_out, samples))
  } else {
    stop(paste0("Unknown pipeline: ", pipeline))
  }
}

#' Transcript quantification
#' @description Calculate the transcript count matrix by parsing the re-alignment file.
#' @param annotation The file path to the annotation file in GFF3 format
#' @param outdir The path to directory to store all output files.
#' @param config Parsed FLAMES configurations.
#' @param pipeline The pipeline type as a character string, either \code{sc_single_sample} (single-cell, single-sample),
#' @param ... Supply sample names as character vector (e.g. \code{samples = c("name1", "name2", ...)}) for muti-sample or bulk pipeline.
#' \code{bulk} (bulk, single or multi-sample), or \code{sc_multi_sample} (single-cell, multiple samples)
#' @return The count matrix will be saved in the output folder as \code{transcript_count.csv.gz}.
#' @examples
#' temp_path <- tempfile()
#' bfc <- BiocFileCache::BiocFileCache(temp_path, ask = FALSE)
#' file_url <- "https://raw.githubusercontent.com/OliverVoogd/FLAMESData/master/data"
#' fastq1 <- bfc[[names(BiocFileCache::bfcadd(bfc, "Fastq1", paste(file_url, "fastq/sample1.fastq.gz", sep = "/")))]]
#' genome_fa <- bfc[[names(BiocFileCache::bfcadd(bfc, "genome.fa", paste(file_url, "SIRV_isoforms_multi-fasta_170612a.fasta", sep = "/")))]]
#' annotation <- bfc[[names(BiocFileCache::bfcadd(bfc, "annot.gtf", paste(file_url, "SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf", sep = "/")))]]
#' outdir <- tempfile()
#' dir.create(outdir)
#' fasta <- annotation_to_fasta(annotation, genome_fa, outdir)
#' config <- jsonlite::fromJSON(create_config(outdir, bambu_isoform_identification = TRUE, min_tr_coverage = 0.1, min_read_coverage = 0.1, min_sup_cnt = 1))
#' file.copy(annotation, file.path(outdir, "isoform_annotated.gtf"))
#' \dontrun{
#' if (!any(is.na(find_bin(c("minimap2", "k8"))))) {
#'   minimap2_realign(
#'     config = config, outdir = outdir,
#'     fq_in = fastq1
#'   )
#'   quantify_transcript_flames(annotation, outdir, config, pipeline = "bulk")
#' }
#' }
#' @export
quantify_transcript <- function(annotation, outdir, config, pipeline = "sc_single_sample", ...) {
  if (config$pipeline_parameters$oarfish_quantification) {
    return(quantify_transcript_oarfish(outdir, config, pipeline, ...))
  } else {
    return(quantify_transcript_flames(annotation, outdir, config, pipeline, ...))
  }
}

# example for Rsamtools
#' @importFrom Rsamtools BamFile  scanBam isIncomplete
# quantify_tmp <- function(bamFileName) {
#    bf <- Rsamtools::BamFile(bamFileName, yieldSize=100)
#    while (Rsamtools::isIncomplete(bf)) {
#        print(Rsamtools::scanBam(bf)[[1]][[1]])
#        Sys.sleep(3)
#    }
# }
