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
#' @return invisible()
#' @import zlibbioc
#' @md
#' @export
find_barcode <- function(fastq, barcodes_file, max_bc_editdistance = 2, max_flank_editdistance = 8,
  reads_out, stats_out, threads = 1, pattern = c(primer = "CTACACGACGCTCTTCCGATCT",
    polyT = paste0(rep("T", 9), collapse = ""), umi_seq = paste0(rep("?", 12),
      collapse = ""), barcode_seq = paste0(rep("?", 16), collapse = "")), full_length_only = FALSE) {
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

  if (sum(grepl("^[35]'TSO$", names(pattern))) == 1) {
    untrimmed_reads <- file.path(
      dirname(reads_out),
      paste0("untrimmed_", basename(reads_out))
    )
    noTSO_reads <- file.path(
      dirname(reads_out),
      paste0("noTSO_", basename(reads_out))
    )
    stopifnot(file.rename(reads_out, untrimmed_reads))
    #cutadapt -a 'TSO' -o reads_out --untrimmed-output noTSO_out in_fq(untrimmed.fq)
    cutadapt(
      c(
        ifelse("3'TSO" %in% names(pattern), "-a", "-g"),
        pattern[[which(grepl("^[35]'TSO$", names(pattern)))]],
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
