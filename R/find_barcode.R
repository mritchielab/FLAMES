#' Match Cell Barcodes
#'
#' @description demultiplex reads with flexiplex
#'
#' @importFrom readr read_delim col_character col_integer col_logical
#' @importFrom dplyr bind_rows
#' @param fastq input FASTQ file path
#' @param barcodes_file path to file containing barcode allow-list, with one barcode in each line
#' @param max_bc_editdistance max edit distances for the barcode sequence
#' @param max_flank_editdistance max edit distances for the flanking sequences (primer and polyT)
#' @param reads_out path of output FASTQ file; if multiple samples are processed,
#' the file name will be appended with the FASTQ file name
#' @param stats_out path of output stats file; if multiple samples are processed,
#' the file name will be appended with the FASTQ file name
#' @param threads number of threads to be used
#' @param full_length_only boolean, when TSO sequence is provided, whether reads without TSO
#' are to be discarded
#' @param pattern named character vector defining the barcode pattern
#' @param TSO_seq TSO sequence to be trimmed
#' @param TSO_prime either 3 (when \code{TSO_seq} is on 3' the end) or 5 (on 5' end)
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' R.utils::gunzip(
#'   filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
#'   destname = bc_allow, remove = FALSE
#' )
#' # single sample
#' find_barcode(
#'   fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
#'   stats_out = file.path(outdir, "bc_stat"),
#'   reads_out = file.path(outdir, "demultiplexed.fq"),
#'   barcodes_file = bc_allow
#' )
#' # multi-sample
#' fastq_dir <- tempfile()
#' dir.create(fastq_dir)
#' file.copy(system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
#'   file.path(fastq_dir, "musc_rps24.fastq.gz"))
#' sampled_lines <- readLines(file.path(fastq_dir, "musc_rps24.fastq.gz"), n = 400)
#' writeLines(sampled_lines, file.path(fastq_dir, "copy.fastq"))
#' result <- find_barcode(
#'   fastq = fastq_dir,
#'   stats_out = file.path(outdir, "bc_stat"),
#'   reads_out = file.path(outdir, "demultiplexed.fq"),
#'   barcodes_file = bc_allow, TSO_seq = "CCCATGTACTCTGCGTTGATACCACTGCTT"
#' )
#' @return a list containing: \code{reads_tb} (tibble of read demultiplexed information) and
#'  \code{input}, \code{output}, \code{read1_with_adapter} from cutadapt report
#'  (if TSO trimming is performed)
#' @md
#' @export
find_barcode <- function(
    fastq, barcodes_file, max_bc_editdistance = 2, max_flank_editdistance = 8,
    reads_out, stats_out, threads = 1, pattern = c(
      primer = "CTACACGACGCTCTTCCGATCT",
      BC = paste0(rep("N", 16), collapse = ""),
      UMI = paste0(rep("N", 12), collapse = ""),
      polyT = paste0(rep("T", 9), collapse = "")
    ), TSO_seq = "", TSO_prime = 3, full_length_only = FALSE) {

  # stats file columns
  # readr's guessing will fail if the file is empty (no reads)
  stats_col_types <- list(
    "Read" = readr::col_character(),
    "CellBarcode" = readr::col_character(),
    "FlankEditDist" = readr::col_integer(),
    "BarcodeEditDist" = readr::col_integer(),
    "UMI" = readr::col_character(),
    "TooShort" = readr::col_logical()
  )

  if (file_test("-f", fastq)) {
    # Single sample
    flexiplex(
      reads_in = fastq, barcodes_file = barcodes_file, bc_as_readid = TRUE,
      max_bc_editdistance = max_bc_editdistance, max_flank_editdistance = max_flank_editdistance,
      pattern = pattern, reads_out = reads_out, stats_out = stats_out, n_threads = threads,
      bc_out = tempfile()
    )
    stats_tb <- readr::read_delim(stats_out, col_types = stats_col_types)
    # for compatibility with multi-sample and TSO trimming
    stats_tb$Outfile <- reads_out
    stats_tb$Sample <- basename(stats_out)
  } else if (file_test("-d", fastq)) {
    # Multi-sample / files
    stats_list <- sapply(list.files(fastq, "\\.(fastq)|(fq)|(fasta)|(fa)$", full.names = TRUE), function(x) {
      stats_out_i <- paste0(stats_out, "_", gsub("\\.gz$", "", basename(x)))
      stats_out_i <- gsub("(\\.fastq$)|(\\.fq$)", "", stats_out_i)
      reads_out_i <- paste0(reads_out, "_", gsub("\\.gz$", "", basename(x)))
      flexiplex(
        reads_in = x, barcodes_file = barcodes_file, bc_as_readid = TRUE,
        max_bc_editdistance = max_bc_editdistance, max_flank_editdistance = max_flank_editdistance,
        pattern = pattern, reads_out = reads_out_i, stats_out = stats_out,
        n_threads = threads, bc_out = tempfile()
      )
      stats_i <- readr::read_delim(stats_out, col_types = stats_col_types)
      stats_i$Outfile <- reads_out_i
      return(stats_i)
    }, simplify = FALSE)
    names(stats_list) <- basename(names(stats_list))
    names(stats_list) <- gsub("\\.(gz)$", "", names(stats_list))
    names(stats_list) <- gsub("(\\.fastq$)|(\\.fq$)", "", names(stats_list))
    stats_tb <- dplyr::bind_rows(stats_list, .id = "Sample") |>
      dplyr::mutate(Sample = factor(Sample), Outfile = factor(Outfile))
  } else {
    stop("The specified path does not exist: ", fastq)
  }

  # Cutadapt TSO trimming
  reports <- list()
  if (is.character(TSO_seq) && nchar(TSO_seq) > 0 && TSO_prime %in% c(3, 5)) {
    # move the original reads to untrimmed_reads
    # run cutadapt and output to reads_out
    # move the untrimmed reads to noTSO_reads if full_length_only
    for (reads_out_i in unique(stats_tb$Outfile)) {
      untrimmed_reads <- file.path(
        dirname(reads_out_i),
        paste0("untrimmed_", basename(reads_out_i))
      )
      noTSO_reads <- file.path(
        dirname(reads_out_i),
        paste0("noTSO_", basename(reads_out_i))
      )
      stopifnot(file.rename(reads_out_i, untrimmed_reads))
      tmp_json_file <- tempfile(fileext = ".json")
      # cutadapt -a 'TSO' -o reads_out --untrimmed-output noTSO_out in_fq(untrimmed.fq)
      cutadapt(
        c(
          ifelse(TSO_prime == 3, "-a", "-g"),
          TSO_seq,
          "-o",
          reads_out_i,
          untrimmed_reads,
          "--json", tmp_json_file,
          if (full_length_only) {
            c("--untrimmed-output", noTSO_reads)
          }
        )
      )
      # parse cutadapt json report
      report <- jsonlite::fromJSON(tmp_json_file)$read_counts
      report <- report[c('read1_with_adapter', 'input', 'output')]
      report$reads_tb <- stats_tb[stats_tb$Outfile == reads_out_i, ]
      reports[[reads_out_i]] <- report
      unlink(tmp_json_file)
    }
  } else {
    cat("Skipping TSO trimming...\n")
    # keep the same output format
    for (reads_out_i in unique(stats_tb$Outfile)) {
      reports[[reads_out_i]] <- list(reads_tb = stats_tb[stats_tb$Outfile == reads_out_i, ])
    }
  }
  if (length(reports) == 1) {
    names(reports) <- basename(names(reports))
  } else {
    names(reports) <- gsub(paste0("^", reads_out, "_?"), "", names(reports))
  }
  return(reports)
}

#' @importFrom utils read.delim
convert_cellranger_bc <- function(bc_allow, bc_from, bc_to) {
  from <- read.delim(bc_from, header = FALSE)$V1
  to <- read.delim(bc_to, header = FALSE)$V1
  allowed <- read.delim(bc_allow, header = FALSE)$V1
  stopifnot(length(from) == length(to))
  stopifnot(all(allowed %in% from))
  return(to[from %in% allowed])
}

#' Plot Cell Barcode demultiplex statistics
#'
#' @description produce a barplot of cell barcode demultiplex statistics
#'
#' @importFrom dplyr bind_rows group_by summarise ungroup count n_distinct
#' @importFrom ggplot2 ggplot aes geom_line scale_y_log10 ylab xlab ggtitle theme
#' @param find_barcode_result output from \code{\link{find_barcode}}
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' fastq_dir <- tempfile()
#' dir.create(fastq_dir)
#' file.copy(system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
#'   file.path(fastq_dir, "musc_rps24.fastq.gz"))
#' sampled_lines <- readLines(file.path(fastq_dir, "musc_rps24.fastq.gz"), n = 400)
#' writeLines(sampled_lines, file.path(fastq_dir, "copy.fastq"))
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' R.utils::gunzip(
#'   filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
#'   destname = bc_allow, remove = FALSE
#' )
#' find_barcode(
#'   fastq = fastq_dir,
#'   stats_out = file.path(outdir, "bc_stat"),
#'   reads_out = file.path(outdir, "demultiplexed.fq"),
#'   barcodes_file = bc_allow, TSO_seq = "CCCATGTACTCTGCGTTGATACCACTGCTT"
#' ) |>
#'   plot_demultiplex()
#' @return a list of ggplot objects:
#' \itemize{
#' \item knee_plot: knee plot of UMI counts before TSO trimming
#' \item flank_editdistance_plot: flanking sequence (adaptor) edit-distance plot
#' \item barcode_editdistance_plot: barcode edit-distance plot
#' \item cutadapt_plot: if TSO trimming is performed, number of reads kept by cutadapt
#' }
#' @md
#' @export
plot_demultiplex <- function(find_barcode_result) {

  knee_plot <- sapply(find_barcode_result, function(x) {
    x$reads_tb
  }, simplify = FALSE) |>
    dplyr::bind_rows() |>
    dplyr::group_by(CellBarcode, Sample) |>
    dplyr::summarise(
      UMI_count = dplyr::n_distinct(UMI)
    ) |>
    dplyr::ungroup() |>
    dplyr::count(UMI_count, Sample) |>
    ggplot2::ggplot(ggplot2::aes(x = UMI_count, y = n, col = Sample)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_log10() +
    ggplot2::ylab("number of cells") +
    ggplot2::xlab("UMI count (before TSO trimming)")
  # remove color legend if only one sample
  if (length(find_barcode_result) == 1) {
    knee_plot <- knee_plot + ggplot2::theme(legend.position = "none")
  }

  flank_editdistance_plot <- sapply(find_barcode_result, function(x) {
    x$reads_tb
  }, simplify = FALSE) |>
    dplyr::bind_rows() |>
    dplyr::group_by(FlankEditDist, Sample) |>
    dplyr::summarise(n_reads = dplyr::n()) |>
    dplyr::ungroup() |>
    ggplot2::ggplot(ggplot2::aes(x = FlankEditDist, y = n_reads, col = Sample)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_log10() +
    ggplot2::ylab("number of reads") +
    ggplot2::xlab("Flanking sequence (adaptor) editdistance")
  if (length(find_barcode_result) == 1) {
    flank_editdistance_plot <- flank_editdistance_plot +
      ggplot2::theme(legend.position = "none")
  }

  # grouped (by sample) barplot
  barcode_editdistance_plot <- sapply(find_barcode_result, function(x) {
    x$reads_tb
  }, simplify = FALSE) |>
    dplyr::bind_rows() |>
    dplyr::group_by(BarcodeEditDist, Sample) |>
    dplyr::ungroup() |>
    ggplot2::ggplot(ggplot2::aes(x = BarcodeEditDist, fill = Sample)) +
    ggplot2::geom_bar(stat = "count", position = "dodge") +
    ggplot2::scale_y_log10() +
    ggplot2::ylab("number of reads") +
    ggplot2::xlab("Barcode editdistance")
  if (length(find_barcode_result) == 1) {
    barcode_editdistance_plot <- barcode_editdistance_plot +
      ggplot2::theme(legend.position = "none")
  }

  # Cutadapt report
  has_cutadapt <- sapply(find_barcode_result, function(x) {
    all(
      "input" %in% names(x),
      "output" %in% names(x),
      "read1_with_adapter" %in% names(x)
    )
  }) |> all()
  if (has_cutadapt) {
    tb <- sapply(find_barcode_result, function(x) {
      c(x$input, x$output, x$read1_with_adapter)
    })
    colnames(tb) <- basename(colnames(tb))
    rownames(tb) <- c("input", "output", "read1_with_adapter")
    tb <- as_tibble(t(tb), rownames = "Sample")
    cutadapt_plot <- tb |>
      tidyr::pivot_longer(
        cols = c(input, output, read1_with_adapter),
        names_to = "Type", values_to = "Count"
      ) |>
      ggplot2::ggplot(ggplot2::aes(x = Sample, y = Count, fill = Type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::ylab("number of reads") +
      ggplot2::xlab("Sample") +
      ggplot2::ggtitle("Cutadapt results")
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  return(list(
    knee_plot = knee_plot,
    flank_editdistance_plot = flank_editdistance_plot,
    barcode_editdistance_plot = barcode_editdistance_plot,
    cutadapt_plot = switch(has_cutadapt, cutadapt_plot, NULL)
  ))
}
