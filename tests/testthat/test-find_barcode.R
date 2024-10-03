test_that("barcode_output_file_identical", {
  outdir <- tempfile()
  dir.create(outdir)
  bc_allow <- file.path(outdir, "bc_allow.tsv")
  R.utils::gunzip(
    filename = system.file("extdata/bc_allow.tsv.gz", package = "FLAMES"),
    destname = bc_allow, remove = FALSE, overwrite = TRUE
  )

  find_barcode(
    max_bc_editdistance = 2, max_flank_editdistance = 8,
    fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
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
    read.delim(file.path(outdir, "stats.tsv"))[, -1]
  )
  expect_identical(
    readLines(system.file('extdata', 'fastq', 'demultiplexed.fq.gz', package = 'FLAMES')),
    readLines(file.path(outdir, "out.fq"), n = 40)
  )
})
test_that("multiple fastq files as one sample", {
  outdirx <- tempfile()
  outdiry <- tempfile()
  fastq_dir <- tempfile()
  c(outdirx, outdiry, fastq_dir) |>
    sapply(dir.create)

  bc_allow <- file.path(outdirx, "bc_allow.tsv")
  R.utils::gunzip(
    filename = system.file("extdata/bc_allow.tsv.gz", package = "FLAMES"),
    destname = bc_allow, remove = FALSE, overwrite = TRUE
  )
  # split the fastq file
  lines <- readLines(system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"))
  i <- 1
  while (length(lines) > 0) {
    writeLines(lines[1:400], file.path(fastq_dir, paste0("musc_rps24_", i, ".fastq")))
    lines <- lines[-(1:400)]
    i <- i + 1
  }

  x <- find_barcode(
    max_bc_editdistance = 2, max_flank_editdistance = 8,
    fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
    barcodes_file = bc_allow,
    reads_out = file.path(outdirx, "out.fq"),
    stats_out = file.path(outdirx, "stats.tsv"),
    threads = 1, pattern = c(
      primer = "CTACACGACGCTCTTCCGATCT",
      BC = paste0(rep("N", 16), collapse = ""),
      UMI = paste0(rep("N", 12), collapse = ""),
      polyT = paste0(rep("T", 9), collapse = "")
    ), TSO_seq = "", TSO_prime = 3, full_length_only = FALSE
  )

  y <- find_barcode(
    max_bc_editdistance = 2, max_flank_editdistance = 8,
    fastq = fastq_dir,
    barcodes_file = bc_allow,
    reads_out = file.path(outdiry, "out.fq"),
    stats_out = file.path(outdiry, "stats.tsv"),
    threads = 1, pattern = c(
      primer = "CTACACGACGCTCTTCCGATCT",
      BC = paste0(rep("N", 16), collapse = ""),
      UMI = paste0(rep("N", 12), collapse = ""),
      polyT = paste0(rep("T", 9), collapse = "")
    ), TSO_seq = "", TSO_prime = 3, full_length_only = FALSE
  )

  expect_identical(
    dplyr::select(x$musc_rps24$reads_tb, -Outfile, -Sample),
    dplyr::select(y[[1]]$reads_tb, -Outfile, -Sample)
  )
  expect_identical(
    read.delim(file.path(outdirx, "stats.tsv")),
    read.delim(file.path(outdiry, "stats.tsv"))
  )
  expect_identical(
    readLines(file.path(outdirx, "out.fq")),
    readLines(file.path(outdiry, "out.fq"))
  )
})

test_that("multiple samples, either file or folder as input", {
  outdirx <- tempfile()
  outdiry <- tempfile()
  fastq_dir <- tempfile()
  c(outdirx, outdiry, fastq_dir) |>
    sapply(dir.create)

  bc_allow <- file.path(outdirx, "bc_allow.tsv")
  R.utils::gunzip(
    filename = system.file("extdata/bc_allow.tsv.gz", package = "FLAMES"),
    destname = bc_allow, remove = FALSE, overwrite = TRUE
  )
  # split the fastq file
  lines <- readLines(system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"))
  i <- 1
  while (length(lines) > 0) {
    writeLines(lines[1:400], file.path(fastq_dir, paste0("musc_rps24_", i, ".fastq")))
    lines <- lines[-(1:400)]
    i <- i + 1
  }

  x <- find_barcode(
    max_bc_editdistance = 2, max_flank_editdistance = 8,
    fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
    barcodes_file = bc_allow,
    reads_out = file.path(outdirx, "out.fq"),
    stats_out = file.path(outdirx, "stats.tsv"),
    threads = 1, pattern = c(
      primer = "CTACACGACGCTCTTCCGATCT",
      BC = paste0(rep("N", 16), collapse = ""),
      UMI = paste0(rep("N", 12), collapse = ""),
      polyT = paste0(rep("T", 9), collapse = "")
    ), TSO_seq = "", TSO_prime = 3, full_length_only = FALSE
  )

  y <- find_barcode(
    max_bc_editdistance = 2, max_flank_editdistance = 8,
    fastq = c(
      'sampleA' = fastq_dir,
      'sampleB' = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
    ),
    barcodes_file = bc_allow,
    reads_out = file.path(outdiry, "out"),
    stats_out = file.path(outdiry, "stats.tsv"),
    threads = 1, pattern = c(
      primer = "CTACACGACGCTCTTCCGATCT",
      BC = paste0(rep("N", 16), collapse = ""),
      UMI = paste0(rep("N", 12), collapse = ""),
      polyT = paste0(rep("T", 9), collapse = "")
    ), TSO_seq = "", TSO_prime = 3, full_length_only = FALSE
  )

  expect_identical(
    dplyr::select(y$sampleA$reads_tb, -Outfile, -Sample),
    dplyr::select(y$sampleB$reads_tb, -Outfile, -Sample)
  )
  expect_identical(
    dplyr::select(y[[1]]$reads_tb, -Outfile, -Sample),
    dplyr::select(x[[1]]$reads_tb, -Outfile, -Sample)
  )
  expect_identical(
    readLines(as.character(y[[1]]$reads_tb$Outfile[1])),
    readLines(as.character(y[[2]]$reads_tb$Outfile[1]))
  )
  expect_identical(
    readLines(as.character(y[[1]]$reads_tb$Outfile[1])),
    readLines(as.character(x[[1]]$reads_tb$Outfile[1]))
  )
})
