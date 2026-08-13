# Pipeline for bulk long read RNA-seq data processing (deprecated)

This function is deprecated. Use
[`BulkPipeline`](https://mritchielab.github.io/FLAMES/reference/BulkPipeline.md)
instead.

## Usage

``` r
bulk_long_pipeline(
  annotation,
  fastq,
  outdir,
  genome_fa,
  minimap2 = NULL,
  config_file
)
```

## Arguments

- annotation:

  The file path to the annotation file in GFF3 / GTF format.

- fastq:

  Path to the FASTQ file or a directory containing FASTQ files. Each
  file will be processed as an individual sample.

- outdir:

  Path to the output directory. If it does not exist, it will be
  created.

- genome_fa:

  The file path to the reference genome in FASTA format.

- minimap2:

  (optional) The path to the minimap2 binary. If not provided, FLAMES
  will use a copy from bioconda via `basilisk`.

- config_file:

  Path to the JSON configuration file. See
  [`create_config`](https://mritchielab.github.io/FLAMES/reference/create_config.md)
  for creating one.

## Value

A `SummarizedExperiment` object containing the transcript counts.

## See also

[`BulkPipeline`](https://mritchielab.github.io/FLAMES/reference/BulkPipeline.md)
for the new pipeline function.
[`SingleCellPipeline`](https://mritchielab.github.io/FLAMES/reference/SingleCellPipeline.md)
for single cell pipelines,
[`MultiSampleSCPipeline`](https://mritchielab.github.io/FLAMES/reference/MultiSampleSCPipeline.md)
for multi sample single cell pipelines.

## Examples

``` r
outdir <- tempfile()
dir.create(outdir)
# simulate 3 samples via sampling
reads <- ShortRead::readFastq(
  system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
)
dir.create(file.path(outdir, "fastq"))
ShortRead::writeFastq(reads[1:100],
  file.path(outdir, "fastq/sample1.fq.gz"),
  mode = "w", full = FALSE
)
reads <- reads[-(1:100)]
ShortRead::writeFastq(reads[1:100],
  file.path(outdir, "fastq/sample2.fq.gz"),
  mode = "w", full = FALSE
)
reads <- reads[-(1:100)]
ShortRead::writeFastq(reads,
  file.path(outdir, "fastq/sample3.fq.gz"),
  mode = "w", full = FALSE
)
# prepare the reference genome
genome_fa <- file.path(outdir, "rps24.fa")
R.utils::gunzip(
  filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
  destname = genome_fa, remove = FALSE
)
se <- bulk_long_pipeline(
  fastq = file.path(outdir, "fastq"),
  annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  outdir = outdir, genome_fa = genome_fa,
  config_file = create_config(outdir, type = "sc_3end", threads = 1, no_flank = TRUE)
)
#> bulk_long_pipeline() is deprecated. Use BulkPipeline() instead.
#> ℹ Writing configuration to: /tmp/RtmpgehRKJ/filebbc05ef7984e/config_file_48064.json
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> FLAMES version 2.7.1 (unknown source)
#> ── Running step: genome_alignment @ Thu Aug 13 10:09:30 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpgehRKJ/filebbc05ef7984e/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpgehRKJ/filebbc05ef7984e/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpgehRKJ/filebbc05ef7984e/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Thu Aug 13 10:09:30 2026 ─────────────
#> ── Running step: read_realignment @ Thu Aug 13 10:09:31 2026 ───────────────────
#> Realigning sample sample1 -> /tmp/RtmpgehRKJ/filebbc05ef7984e/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/RtmpgehRKJ/filebbc05ef7984e/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/RtmpgehRKJ/filebbc05ef7984e/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Thu Aug 13 10:09:31 2026 ──────────
#> Pipeline saved to /tmp/RtmpgehRKJ/filebbc05ef7984e/pipeline.rds
se
#> class: SummarizedExperiment 
#> dim: 10 3 
#> metadata(0):
#> assays(1): counts
#> rownames(10): ENSMUSG00000025290.17_19_5159_1
#>   ENSMUSG00000025290.17_19_5159_2 ... ENSMUST00000169826.2
#>   ENSMUST00000225023.1
#> rowData names(0):
#> colnames(3): sample1 sample2 sample3
#> colData names(0):
```
