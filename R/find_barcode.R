#' Match Cell Barcodes
#'
#' @description demultiplex reads with flexiplex
#'
#' @param fastq input FASTQ file path
#' @param barcodes_file path to file containing barcode allow-list, with one barcode in each line
#' @param max_bc_editdistance max edit distances for the barcode sequence
#' @param max_flank_editdistance max edit distances for the flanking sequences (primer and polyT)
#' @param reads_out path of output FASTQ file
#' @param stats_out path of output stats file
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
#' R.utils::gunzip(filename = system.file("extdata/bc_allow.tsv.gz", package = "FLAMES"), destname = bc_allow, remove = FALSE)
#' find_barcode(
#'   fastq = system.file("extdata/fastq", package = "FLAMES"),
#'   stats_out = file.path(outdir, "bc_stat"),
#'   reads_out = file.path(outdir, "demultiplexed.fq.gz"),
#'   barcodes_file = bc_allow
#' )
#' @return invisible()
#' @import zlibbioc
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
  if (file_test("-f", fastq)) {
    flexiplex(
      reads_in = fastq, barcodes_file = barcodes_file, bc_as_readid = TRUE,
      max_bc_editdistance = max_bc_editdistance, max_flank_editdistance = max_flank_editdistance,
      pattern = pattern, reads_out = reads_out, stats_out = stats_out, n_threads = threads,
      bc_out = tempfile()
    )
  } else if (file_test("-d", fastq)) {
    sapply(list.files(fastq, "\\.(fastq)|(fq)|(fasta)|(fa)$", full.names = TRUE), function(x) {
      flexiplex(
        reads_in = x, barcodes_file = barcodes_file, bc_as_readid = TRUE,
        max_bc_editdistance = max_bc_editdistance, max_flank_editdistance = max_flank_editdistance,
        pattern = pattern, reads_out = reads_out, stats_out = stats_out,
        n_threads = threads, bc_out = tempfile()
      )
    })
  } else {
    stop("The specified path does not exist: ", fastq)
  }

  if (is.character(TSO_seq) && nchar(TSO_seq) > 0 && TSO_prime %in% c(3, 5)) {
    untrimmed_reads <- file.path(
      dirname(reads_out),
      paste0("untrimmed_", basename(reads_out))
    )
    noTSO_reads <- file.path(
      dirname(reads_out),
      paste0("noTSO_", basename(reads_out))
    )
    stopifnot(file.rename(reads_out, untrimmed_reads))
    # cutadapt -a 'TSO' -o reads_out --untrimmed-output noTSO_out in_fq(untrimmed.fq)
    cutadapt(
      c(
        ifelse(TSO_prime == 3, "-a", "-g"),
        TSO_seq,
        "-o",
        reads_out,
        untrimmed_reads,
        if (full_length_only) {
          c("--untrimmed-output", noTSO_reads)
        }
      )
    )
  } else {
    cat("Skipping TSO trimming...\n")
  }
}

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
#' @param outdir folder containing the \code{matched_barcode_stat} file, or
#' \code{matched_barcode_stat.SAMPLE} files. Ignored if \code{stats_file} is provided.
#' @param stats_file \code{matched_barcode_stat} file(s) from which the statistics to be plotted.
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' bc_allow <- file.path(outdir, "bc_allow.tsv")
#' R.utils::gunzip(filename = system.file("extdata/bc_allow.tsv.gz", package = "FLAMES"), destname = bc_allow, remove = FALSE)
#' find_barcode(
#'   fastq = system.file("extdata/fastq", package = "FLAMES"),
#'   stats_out = file.path(outdir, "bc_stat"),
#'   reads_out = file.path(outdir, "demultiplexed.fq.gz"),
#'   barcodes_file = bc_allow
#' )
#' plot_demultiplex(stats_file = file.path(outdir, "bc_stat"))
#' @return a \code{ggplot} object of the barcode plot
#' @md
#' @export
plot_demultiplex <- function(outdir, stats_file) {
  if (missing(stats_file)) {
    stats_file <- list.files(path = outdir, pattern = "matched_barcode_stat.*", full.names = TRUE)
  }
  stopifnot("matched_barcode_stat file not found under outdir! Please specifiy them with 'stats_file' arguement" = length(stats_file) >= 1)

  if (length(stats_file) == 1) {
    stats_df <- read.delim(stats_file)
  } else { # handle multi-samples: append file name to df$CellBarcode
    stats_dfs <- sapply(stats_file, read.delim, simplify = FALSE)
    for (sample in names(stats_dfs)) {
      stats_dfs[[sample]]$CellBarcode <- paste(stats_dfs[[sample]]$CellBarcode, sample, sep = "-")
    }
    stats_df <- do.call(rbind, stats_dfs)
  }

  knee_plot <- stats_df$CellBarcode |>
    table() |>
    data.frame() |>
    dplyr::arrange(dplyr::desc(Freq)) %>%
    cbind(x = seq_len(nrow(.))) |>
    ggplot2::ggplot(ggplot2::aes(x = x, y = Freq)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_log10()

  editdistance_plot <- stats_df[, c("FlankEditDist", "BarcodeEditDist")] |>
    table() |>
    tidyr::as_tibble() |>
    ggplot2::ggplot(ggplot2::aes(x = FlankEditDist, y = n, fill = BarcodeEditDist)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::ylab("number of reads") +
    ggplot2::xlab("Adaptor editdistance")

  return(list(
    knee_plot = knee_plot,
    editdistance_plot = editdistance_plot
  ))
}
