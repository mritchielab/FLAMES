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
#' @param pattern named character vector defining the barcode pattern
#' @examples
#' outdir <- tempfile()
#' dir.create(outdir)
#' bc_allow <- file.path(outdir, 'bc_allow.tsv')
#' R.utils::gunzip(filename = system.file('extdata/bc_allow.tsv.gz', package = 'FLAMES'), destname = bc_allow, remove = FALSE)
#' find_barcode(
#'    fastq = system.file('extdata/fastq', package = 'FLAMES'),
#'    stats_out = file.path(outdir, 'bc_stat'),
#'    reads_out = file.path(outdir, 'demultiplexed.fq.gz'),
#'    barcodes_file = bc_allow
#')
#' @md
#' @export
find_barcode <- function(fastq, barcodes_file, max_bc_editdistance = 2, max_flank_editdistance = 8,
  reads_out, stats_out, threads = 1, pattern = c(primer = "CTACACGACGCTCTTCCGATCT",
    polyT = paste0(rep("T", 9), collapse = ""), umi_seq = paste0(rep("?", 12),
      collapse = ""), barcode_seq = paste0(rep("?", 16), collapse = ""))) {
  if (file_test("-f", fastq)) {
    flexiplex(reads_in = fastq, barcodes_file = barcodes_file, bc_as_readid = TRUE,
      max_bc_editdistance = max_bc_editdistance, max_flank_editdistance = max_flank_editdistance,
      pattern = pattern, reads_out = reads_out, stats_out = stats_out, n_threads = threads,
      bc_out = tempfile())
  } else if (file_test("-d", fastq)) {
    sapply(list.files(fastq, "\\.(fastq)|(fq)|(fasta)|(fa)$", full.names=TRUE), function(x) {
      flexiplex(reads_in = x, barcodes_file = barcodes_file, bc_as_readid = TRUE,
        max_bc_editdistance = max_bc_editdistance, max_flank_editdistance = max_flank_editdistance,
        pattern = pattern, reads_out = reads_out, stats_out = stats_out,
        n_threads = threads, bc_out = tempfile())
    })
  } else {
    stop("The specified path does not exist: ", fastq)
  }
}
