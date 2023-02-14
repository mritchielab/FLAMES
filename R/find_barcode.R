#' Match Cell Barcodes
#'
#' @description Match cell barcodes in the given fastq directory with the given barcode allow-list. For each read, the left flanking
#' sequence is located, FLAMES then takes the next 16 characters and match it to barcodes in the allow-list. If there is an unambigious match
#' within the given edit distance (\code{MAX_DIST}), the barcode and following \code{UMI_LEN} characters are tirmmed, along with potential
#' polyT tail. The trimmed read is then saved to \code{out_fastq}, with the identier field formatted as \[barcode\]_\[UMI\]#\[original read ID\].
#' 
#' @param fastq_dir directory containing fastq files to match
#' @param stats_file NEEDED
#' @param out_fastq output filename for matched barcodes
#' @param ref_csv NEEDED
#' @param MAX_DIST int; maximum edit distance
#' @param UMI_LEN int; length of UMI sequences
#' @param left_seq String; sequence that appears at the left of the barcode
#' @param min_length int; minimum read length to be filtered after timming barcodes
#' @param reverse_complement boolean; whether to check the reverse complement of the reads
#' @param fixed_range boolean; deprecated, whether to skip finding flanking sequence by infering
#' its position from previous reads. Setting to \code{TRUE} may decrease performance and accuracy.
#'
#' @return returns a list containing statistics of the reads demultiplexed.
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' bc_allow <- file.path(outdir, 'bc_allow.tsv')
#' R.utils::gunzip(filename = system.file('extdata/bc_allow.tsv.gz', package = 'FLAMES'), destname = bc_allow, remove = FALSE)
#' find_barcode(
#'    fastq_dir = system.file('extdata/fastq', package = 'FLAMES'),
#'    stats_file = file.path(outdir, 'bc_stat'),
#'    out_fastq = file.path(outdir, 'demultiplexed.fq.gz'),
#'    ref_csv = bc_allow,
#'    MAX_DIST = 2,
#'    UMI_LEN = 10
#')
#' @md
#' @export
find_barcode <- function(fastq_dir, stats_file, out_fastq, ref_csv, MAX_DIST, UMI_LEN = 10L,
  left_seq = "CTACACGACGCTCTTCCGATCT", min_length = 20L, reverse_complement = TRUE,
  fixed_range = FALSE) {
  stats <- match_cell_barcode(fastq_dir, stats_file, out_fastq, ref_csv, MAX_DIST,
    UMI_LEN, left_seq, min_length, reverse_complement, fixed_range)
  demultiplexed_counts <- read.csv(stats_file, row.names = "cell_barcode")
  return(list(demultiplex_stats = as.data.frame(stats), column_data = demultiplexed_counts))
}

#' Barcode demultiplexing QC plots
#' @description Plot the barcode demultiplexing statistics
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr '%>%'
#' @importFrom dplyr filter
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom ggplot2 ggplot aes geom_bar ggtitle coord_polar element_blank theme position_stack theme_bw geom_histogram ggtitle ylab xlab geom_text
#' @param sce The \code{SingleCellExperiment} object from FLAMES pipeline, or the returned list from \code{find_barcode}
#' @return a list of QC plots for the barcode demultiplexing step (\code{find_barcode})
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' bc_allow <- file.path(outdir, 'bc_allow.tsv')
#' R.utils::gunzip(filename = system.file('extdata/bc_allow.tsv.gz', package = 'FLAMES'), destname = bc_allow, remove = FALSE)
#' barcode_info <- find_barcode(
#'    fastq_dir = system.file('extdata/fastq', package = 'FLAMES'),
#'    stats_file = file.path(outdir, 'bc_stat'),
#'    out_fastq = file.path(outdir, 'demultiplexed.fq.gz'),
#'    ref_csv = bc_allow,
#'    MAX_DIST = 2,
#'    UMI_LEN = 10
#')
#' barcode_info_plots(barcode_info)
#' @md
#' @export
barcode_info_plots <- function(sce) {
  if (is.list(sce)) {
    demultiplex_stats <- sce$demultiplex_stats
    demultiplexed_counts <- sce$column_data[, "reads_demultiplexed", drop = F]
  } else {
    demultiplex_stats <- sce@metadata$demultiplex_stats
    demultiplexed_counts <- colData(sce)[, "reads_demultiplexed", drop = F]
  }
  summary_pie <- tidyr::pivot_longer(demultiplex_stats, everything()) %>%
    dplyr::filter(!name %in% c("total")) %>%
    ggplot(aes(x = "", y = value, label = value, fill = name)) + geom_bar(stat = "identity") +
    geom_text(position = position_stack(vjust = 0.5)) + coord_polar("y") + labs(x = NULL,
    y = NULL) + ggtitle("Reads demultiplexed summary") + theme_bw() + theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank(),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank())
  distribution_hist <- ggplot(demultiplexed_counts, aes(x = reads_demultiplexed)) +
    geom_histogram(bins = ifelse(length(table(demultiplexed_counts)) < 20, length(table(demultiplexed_counts)),
      20)) + xlab("Reads demultiplexed") + ylab("Cell barcodes") + ggtitle("Reads demultiplexed distribution")
  return(list(summary_pie = summary_pie, distribution_hist = distribution_hist))
}
