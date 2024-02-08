test_that("barcode_output_file_identical", {
  outdir <- tempdir()
  bc_allow <- file.path(outdir, "bc_allow.tsv")
  R.utils::gunzip(
    filename = system.file("extdata/bc_allow.tsv.gz", package = "FLAMES"),
    destname = bc_allow, remove = FALSE
  )

  find_barcode(
    max_bc_editdistance = 2, max_flank_editdistance = 8,
    fastq = system.file("extdata/fastq", package = "FLAMES"),
    barcodes_file = bc_allow,
    reads_out = file.path(outdir, "out.fq"),
    stats_out = file.path(outdir, "stats.tsv"),
    threads = 1, pattern = c(
      primer = "CTACACGACGCTCTTCCGATCT",
      BC = paste0(rep("N", 16), collapse = ""),
      UMI = paste0(rep("N", 12), collapse = ""),
      polyT = paste0(rep("T", 9), collapse = "")
    ), TSO_seq = "", TSO_prime = 3, full_length_only = FALSE
  )

  expect_identical(
    read.delim(test_path("bc_stat")),
    read.delim(file.path(outdir, "stats.tsv"))
  )
  expect_identical(
    read.delim(test_path("demultiplexed.fq")),
    read.delim(file.path(outdir, "out.fq"))
  )
})
