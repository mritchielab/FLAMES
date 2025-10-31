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
#> Writing configuration parameters to:  /tmp/Rtmpbzssfl/filea709406a8d82/config_file_42761.json 
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> ── Running step: barcode_demultiplex @ Fri Oct 31 06:02:56 2025 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> FLEXIPLEX 0.96.2
#> Setting max barcode edit distance to 2
#> Setting max flanking sequence edit distance to 8
#> Setting read IDs to be  replaced
#> Setting number of threads to 8
#> Search pattern: 
#> primer: CTACACGACGCTCTTCCGATCT
#> BC: NNNNNNNNNNNNNNNN
#> UMI: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Setting known barcodes from /tmp/Rtmpbzssfl/filea709406a8d82/bc_allow.tsv
#> Number of known barcodes: 143
#> Processing file: /tmp/Rtmpbzssfl/filea709406a8d82/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/Rtmpbzssfl/filea709406a8d82/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/Rtmpbzssfl/filea709406a8d82/fastq/sample3.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 393
#> Number of reads where at least one barcode was found: 368
#> Number of reads with exactly one barcode match: 364
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
#> FLEXIPLEX 0.96.2
#> Setting max barcode edit distance to 2
#> Setting max flanking sequence edit distance to 8
#> Setting read IDs to be  replaced
#> Setting number of threads to 8
#> Search pattern: 
#> primer: CTACACGACGCTCTTCCGATCT
#> BC: NNNNNNNNNNNNNNNN
#> UMI: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Setting known barcodes from /tmp/Rtmpbzssfl/filea709406a8d82/bc_allow.tsv
#> Number of known barcodes: 143
#> Processing file: /tmp/Rtmpbzssfl/filea709406a8d82/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 100
#> Number of reads where at least one barcode was found: 92
#> Number of reads with exactly one barcode match: 91
#> Number of chimera reads: 1
#> All done!
#> Reads    Barcodes
#> 4    1
#> 3    9
#> 2    9
#> 1    44
#> FLEXIPLEX 0.96.2
#> Setting max barcode edit distance to 2
#> Setting max flanking sequence edit distance to 8
#> Setting read IDs to be  replaced
#> Setting number of threads to 8
#> Search pattern: 
#> primer: CTACACGACGCTCTTCCGATCT
#> BC: NNNNNNNNNNNNNNNN
#> UMI: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Setting known barcodes from /tmp/Rtmpbzssfl/filea709406a8d82/bc_allow.tsv
#> Number of known barcodes: 143
#> Processing file: /tmp/Rtmpbzssfl/filea709406a8d82/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 100
#> Number of reads where at least one barcode was found: 95
#> Number of reads with exactly one barcode match: 94
#> Number of chimera reads: 0
#> All done!
#> Reads    Barcodes
#> 4    2
#> 3    3
#> 2    16
#> 1    47
#> FLEXIPLEX 0.96.2
#> Setting max barcode edit distance to 2
#> Setting max flanking sequence edit distance to 8
#> Setting read IDs to be  replaced
#> Setting number of threads to 8
#> Search pattern: 
#> primer: CTACACGACGCTCTTCCGATCT
#> BC: NNNNNNNNNNNNNNNN
#> UMI: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Setting known barcodes from /tmp/Rtmpbzssfl/filea709406a8d82/bc_allow.tsv
#> Number of known barcodes: 143
#> Processing file: /tmp/Rtmpbzssfl/filea709406a8d82/fastq/sample3.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 193
#> Number of reads where at least one barcode was found: 181
#> Number of reads with exactly one barcode match: 179
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
#> ── Running step: genome_alignment @ Fri Oct 31 06:02:57 2025 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/Rtmpbzssfl/filea709406a8d82/sampleA_matched_reads.fastq.gz -> /tmp/Rtmpbzssfl/filea709406a8d82/sampleA_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> Aligning sample /tmp/Rtmpbzssfl/filea709406a8d82/sample1_matched_reads.fastq.gz -> /tmp/Rtmpbzssfl/filea709406a8d82/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> Aligning sample /tmp/Rtmpbzssfl/filea709406a8d82/sample2_matched_reads.fastq.gz -> /tmp/Rtmpbzssfl/filea709406a8d82/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> Aligning sample /tmp/Rtmpbzssfl/filea709406a8d82/sample3_matched_reads.fastq.gz -> /tmp/Rtmpbzssfl/filea709406a8d82/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: gene_quantification @ Fri Oct 31 06:02:57 2025 ────────────────
#> 06:02:57 AM Fri Oct 31 2025 quantify genes 
#> Using BAM(s): /tmp/Rtmpbzssfl/filea709406a8d82/sampleA_align2genome.bam,
#> /tmp/Rtmpbzssfl/filea709406a8d82/sample1_align2genome.bam,
#> /tmp/Rtmpbzssfl/filea709406a8d82/sample2_align2genome.bam, and
#> /tmp/Rtmpbzssfl/filea709406a8d82/sample3_align2genome.bam
#> ── Running step: isoform_identification @ Fri Oct 31 06:02:59 2025 ─────────────
#> ── Running step: read_realignment @ Fri Oct 31 06:02:59 2025 ───────────────────
#> Checking for fastq file(s) /tmp/Rtmpbzssfl/filea709406a8d82/fastq, /tmp/Rtmpbzssfl/filea709406a8d82/fastq/sample1.fq.gz, /tmp/Rtmpbzssfl/filea709406a8d82/fastq/sample2.fq.gz, /tmp/Rtmpbzssfl/filea709406a8d82/fastq/sample3.fq.gz
#>  files found
#> Checking for fastq file(s) /tmp/Rtmpbzssfl/filea709406a8d82/sampleA_matched_reads.fastq.gz, /tmp/Rtmpbzssfl/filea709406a8d82/sample1_matched_reads.fastq.gz, /tmp/Rtmpbzssfl/filea709406a8d82/sample2_matched_reads.fastq.gz, /tmp/Rtmpbzssfl/filea709406a8d82/sample3_matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/Rtmpbzssfl/filea709406a8d82/sampleA_matched_reads_dedup.fastq.gz, /tmp/Rtmpbzssfl/filea709406a8d82/sample1_matched_reads_dedup.fastq.gz, /tmp/Rtmpbzssfl/filea709406a8d82/sample2_matched_reads_dedup.fastq.gz, /tmp/Rtmpbzssfl/filea709406a8d82/sample3_matched_reads_dedup.fastq.gz
#>  files found
#> Realigning sample /tmp/Rtmpbzssfl/filea709406a8d82/sampleA_matched_reads_dedup.fastq.gz -> /tmp/Rtmpbzssfl/filea709406a8d82/sampleA_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> Realigning sample /tmp/Rtmpbzssfl/filea709406a8d82/sample1_matched_reads_dedup.fastq.gz -> /tmp/Rtmpbzssfl/filea709406a8d82/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> Realigning sample /tmp/Rtmpbzssfl/filea709406a8d82/sample2_matched_reads_dedup.fastq.gz -> /tmp/Rtmpbzssfl/filea709406a8d82/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> Realigning sample /tmp/Rtmpbzssfl/filea709406a8d82/sample3_matched_reads_dedup.fastq.gz -> /tmp/Rtmpbzssfl/filea709406a8d82/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> ── Running step: transcript_quantification @ Fri Oct 31 06:03:00 2025 ──────────
#> Pipeline saved to /tmp/Rtmpbzssfl/filea709406a8d82/pipeline.rds
```
