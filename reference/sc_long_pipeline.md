# Pipeline for Single Cell Data (deprecated)

This function is deprecated. Please use \[SingleCellPipeline()\]
instead.

## Usage

``` r
sc_long_pipeline(
  annotation,
  fastq,
  outdir,
  genome_fa,
  minimap2 = NULL,
  barcodes_file = NULL,
  expect_cell_number = NULL,
  config_file = NULL
)
```

## Arguments

- annotation:

  The file path to the annotation file in GFF3 format

- fastq:

  The file path to input fastq file

- outdir:

  The path to directory to store all output files.

- genome_fa:

  The file path to genome fasta file.

- minimap2:

  Path to minimap2, optional.

- barcodes_file:

  The file with expected cell barcodes, with each barcode on a new line.

- expect_cell_number:

  The expected number of cells in the sample. This is used if
  `barcodes_file` is not provided. See `BLAZE` for more details.

- config_file:

  File path to the JSON configuration file.

## Value

A `SingleCellPipeline` object containing the transcript counts.

## See also

[`SingleCellPipeline`](https://mritchielab.github.io/FLAMES/reference/SingleCellPipeline.md)
for the new pipeline interface,
[`BulkPipeline`](https://mritchielab.github.io/FLAMES/reference/BulkPipeline.md)
for bulk long data,
[`MultiSampleSCPipeline`](https://mritchielab.github.io/FLAMES/reference/MultiSampleSCPipeline.md)
for multi sample single cell pipelines.

## Examples

``` r
outdir <- tempfile()
dir.create(outdir)
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
sce <- FLAMES::sc_long_pipeline(
  genome_fa = genome_fa,
  fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
  annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  outdir = outdir,
  barcodes_file = bc_allow,
  config_file = FLAMES::create_config(
    outdir,
    pipeline_parameters.demultiplexer = "flexiplex"
  )
)
#> Writing configuration parameters to:  /tmp/RtmpM0tt18/filea709269a191f/config_file_42761.json 
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> ── Running step: barcode_demultiplex @ Fri Oct 31 05:41:39 2025 ────────────────
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
#> Setting known barcodes from /tmp/RtmpM0tt18/filea709269a191f/bc_allow.tsv
#> Number of known barcodes: 143
#> Processing file: /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
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
#> ── Running step: genome_alignment @ Fri Oct 31 05:41:39 2025 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpM0tt18/filea709269a191f/matched_reads.fastq.gz -> /tmp/RtmpM0tt18/filea709269a191f/align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: gene_quantification @ Fri Oct 31 05:41:40 2025 ────────────────
#> 05:41:40 AM Fri Oct 31 2025 quantify genes 
#> Using BAM(s): /tmp/RtmpM0tt18/filea709269a191f/align2genome.bam
#> ── Running step: isoform_identification @ Fri Oct 31 05:41:40 2025 ─────────────
#> ── Running step: read_realignment @ Fri Oct 31 05:41:40 2025 ───────────────────
#> Checking for fastq file(s) /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpM0tt18/filea709269a191f/matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpM0tt18/filea709269a191f/matched_reads_dedup.fastq.gz
#>  files found
#> Realigning sample /tmp/RtmpM0tt18/filea709269a191f/matched_reads_dedup.fastq.gz -> /tmp/RtmpM0tt18/filea709269a191f/realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> ── Running step: transcript_quantification @ Fri Oct 31 05:41:40 2025 ──────────
#> Pipeline saved to /tmp/RtmpM0tt18/filea709269a191f/pipeline.rds
```
