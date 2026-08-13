# Pipeline for Multi-sample Single Cell Data (deprecated)

This function is deprecated. Please use
[`MultiSampleSCPipeline`](https://mritchielab.github.io/FLAMES/reference/MultiSampleSCPipeline.md).

## Usage

``` r
sc_long_multisample_pipeline(
  annotation,
  fastqs,
  outdir,
  genome_fa,
  minimap2 = NULL,
  barcodes_file = NULL,
  expect_cell_numbers = NULL,
  config_file = NULL
)
```

## Arguments

- annotation:

  The file path to the annotation file in GFF3 format

- fastqs:

  The file path to input fastq file

- outdir:

  The path to directory to store all output files.

- genome_fa:

  The file path to genome fasta file.

- minimap2:

  Path to minimap2, optional.

- barcodes_file:

  The file with expected cell barcodes, with each barcode on a new line.

- expect_cell_numbers:

  The expected number of cells in the sample. This is used if
  `barcodes_file` is not provided. See `BLAZE` for more details.

- config_file:

  File path to the JSON configuration file.

## Value

A list of `SingleCellExperiment` objects, one for each sample.

## See also

[`MultiSampleSCPipeline`](https://mritchielab.github.io/FLAMES/reference/MultiSampleSCPipeline.md)
for the new pipeline interface,
[`SingleCellPipeline`](https://mritchielab.github.io/FLAMES/reference/SingleCellPipeline.md)
for single-sample pipeline,
[`BulkPipeline`](https://mritchielab.github.io/FLAMES/reference/BulkPipeline.md)
for bulk long data.

## Examples

``` r
reads <- ShortRead::readFastq(
  system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
)
outdir <- tempfile()
dir.create(outdir)
dir.create(file.path(outdir, "fastq"))
bc_allow <- file.path(outdir, "bc_allow.tsv")
genome_fa <- file.path(outdir, "rps24.fa")
R.utils::gunzip(
  filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
  destname = bc_allow, remove = FALSE
)
R.utils::gunzip(
  filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
  destname = genome_fa, remove = FALSE
)
ShortRead::writeFastq(reads[1:100],
  file.path(outdir, "fastq/sample1.fq.gz"), mode = "w", full = FALSE)
reads <- reads[-(1:100)]
ShortRead::writeFastq(reads[1:100],
  file.path(outdir, "fastq/sample2.fq.gz"), mode = "w", full = FALSE)
reads <- reads[-(1:100)]
ShortRead::writeFastq(reads,
  file.path(outdir, "fastq/sample3.fq.gz"), mode = "w", full = FALSE)

sce_list <- FLAMES::sc_long_multisample_pipeline(
  annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  fastqs = c("sampleA" = file.path(outdir, "fastq"),
    "sample1" = file.path(outdir, "fastq", "sample1.fq.gz"),
    "sample2" = file.path(outdir, "fastq", "sample2.fq.gz"),
    "sample3" = file.path(outdir, "fastq", "sample3.fq.gz")),
  outdir = outdir,
  genome_fa = genome_fa,
  barcodes_file = rep(bc_allow, 4),
  config_file = create_config(
    outdir,
    pipeline_parameters.demultiplexer = "flexiplex"
  )
)
#> sc_long_multisample_pipeline is deprecated, please use MultiSampleSCPipeline instead.
#> ℹ Writing configuration to: /tmp/RtmpgehRKJ/filebbc05d5e422f/config_file_48064.json
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> FLAMES version 2.7.1 (unknown source)
#> ── Running step: barcode_demultiplex @ Thu Aug 13 10:10:59 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/RtmpgehRKJ/filebbc05d5e422f/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 8
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> CB:Z: tag field: CB
#> Processing file: /tmp/RtmpgehRKJ/filebbc05d5e422f/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/RtmpgehRKJ/filebbc05d5e422f/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/RtmpgehRKJ/filebbc05d5e422f/fastq/sample3.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 393
#> Number of reads where at least one barcode was found: 368
#> Number of chimera reads: 1
#> All done!
#> Reads    Barcodes
#> 10   2
#> 9    2
#> 8    5
#> 7    4
#> 6    3
#> 5    7
#> 4    14
#> 3    14
#> 2    29
#> 1    57
#> Loading known barcodes from /tmp/RtmpgehRKJ/filebbc05d5e422f/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 8
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> CB:Z: tag field: CB
#> Processing file: /tmp/RtmpgehRKJ/filebbc05d5e422f/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 100
#> Number of reads where at least one barcode was found: 92
#> Number of chimera reads: 1
#> All done!
#> Reads    Barcodes
#> 4    1
#> 3    9
#> 2    9
#> 1    44
#> Loading known barcodes from /tmp/RtmpgehRKJ/filebbc05d5e422f/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 8
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> CB:Z: tag field: CB
#> Processing file: /tmp/RtmpgehRKJ/filebbc05d5e422f/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 100
#> Number of reads where at least one barcode was found: 95
#> Number of chimera reads: 0
#> All done!
#> Reads    Barcodes
#> 4    2
#> 3    3
#> 2    16
#> 1    47
#> Loading known barcodes from /tmp/RtmpgehRKJ/filebbc05d5e422f/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 8
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> CB:Z: tag field: CB
#> Processing file: /tmp/RtmpgehRKJ/filebbc05d5e422f/fastq/sample3.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 193
#> Number of reads where at least one barcode was found: 181
#> Number of chimera reads: 0
#> All done!
#> Reads    Barcodes
#> 7    1
#> 6    1
#> 5    1
#> 4    7
#> 3    10
#> 2    27
#> 1    53
#> ── Running step: genome_alignment @ Thu Aug 13 10:11:00 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpgehRKJ/filebbc05d5e422f/sampleA_matched_reads.fastq.gz -> /tmp/RtmpgehRKJ/filebbc05d5e422f/sampleA_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> Aligning sample /tmp/RtmpgehRKJ/filebbc05d5e422f/sample1_matched_reads.fastq.gz -> /tmp/RtmpgehRKJ/filebbc05d5e422f/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> Aligning sample /tmp/RtmpgehRKJ/filebbc05d5e422f/sample2_matched_reads.fastq.gz -> /tmp/RtmpgehRKJ/filebbc05d5e422f/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> Aligning sample /tmp/RtmpgehRKJ/filebbc05d5e422f/sample3_matched_reads.fastq.gz -> /tmp/RtmpgehRKJ/filebbc05d5e422f/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: gene_quantification @ Thu Aug 13 10:11:01 2026 ────────────────
#> 10:11:01 AM Thu Aug 13 2026 quantify genes 
#> Using BAM(s): /tmp/RtmpgehRKJ/filebbc05d5e422f/sampleA_align2genome.bam,
#> /tmp/RtmpgehRKJ/filebbc05d5e422f/sample1_align2genome.bam,
#> /tmp/RtmpgehRKJ/filebbc05d5e422f/sample2_align2genome.bam, and
#> /tmp/RtmpgehRKJ/filebbc05d5e422f/sample3_align2genome.bam
#> ── Running step: isoform_identification @ Thu Aug 13 10:11:02 2026 ─────────────
#> ── Running step: read_realignment @ Thu Aug 13 10:11:02 2026 ───────────────────
#> Checking for fastq file(s) /tmp/RtmpgehRKJ/filebbc05d5e422f/fastq, /tmp/RtmpgehRKJ/filebbc05d5e422f/fastq/sample1.fq.gz, /tmp/RtmpgehRKJ/filebbc05d5e422f/fastq/sample2.fq.gz, /tmp/RtmpgehRKJ/filebbc05d5e422f/fastq/sample3.fq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpgehRKJ/filebbc05d5e422f/sampleA_matched_reads.fastq.gz, /tmp/RtmpgehRKJ/filebbc05d5e422f/sample1_matched_reads.fastq.gz, /tmp/RtmpgehRKJ/filebbc05d5e422f/sample2_matched_reads.fastq.gz, /tmp/RtmpgehRKJ/filebbc05d5e422f/sample3_matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpgehRKJ/filebbc05d5e422f/sampleA_matched_reads_dedup.fastq.gz, /tmp/RtmpgehRKJ/filebbc05d5e422f/sample1_matched_reads_dedup.fastq.gz, /tmp/RtmpgehRKJ/filebbc05d5e422f/sample2_matched_reads_dedup.fastq.gz, /tmp/RtmpgehRKJ/filebbc05d5e422f/sample3_matched_reads_dedup.fastq.gz
#>  files found
#> Realigning sample /tmp/RtmpgehRKJ/filebbc05d5e422f/sampleA_matched_reads_dedup.fastq.gz -> /tmp/RtmpgehRKJ/filebbc05d5e422f/sampleA_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> Realigning sample /tmp/RtmpgehRKJ/filebbc05d5e422f/sample1_matched_reads_dedup.fastq.gz -> /tmp/RtmpgehRKJ/filebbc05d5e422f/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> Realigning sample /tmp/RtmpgehRKJ/filebbc05d5e422f/sample2_matched_reads_dedup.fastq.gz -> /tmp/RtmpgehRKJ/filebbc05d5e422f/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> Realigning sample /tmp/RtmpgehRKJ/filebbc05d5e422f/sample3_matched_reads_dedup.fastq.gz -> /tmp/RtmpgehRKJ/filebbc05d5e422f/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> ── Running step: transcript_quantification @ Thu Aug 13 10:11:03 2026 ──────────
#> Pipeline saved to /tmp/RtmpgehRKJ/filebbc05d5e422f/pipeline.rds
```
