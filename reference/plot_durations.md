# Plot pipeline step durations

This function creates a horizontal bar plot showing the duration of each
pipeline step using ggplot2.

## Usage

``` r
plot_durations(x)

# S4 method for class 'FLAMES.Pipeline'
plot_durations(x)
```

## Arguments

- x:

  A FLAMES.Pipeline object.

## Value

A ggplot2 object.

## Examples

``` r
pipeline <- example_pipeline("BulkPipeline")
#> Writing configuration parameters to:  /tmp/RtmpjRxi1E/file95d242b7da9b/config_file_38354.json 
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- run_FLAMES(pipeline)
#> ── Running step: genome_alignment @ Thu Feb 19 01:56:09 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpjRxi1E/file95d242b7da9b/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpjRxi1E/file95d242b7da9b/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpjRxi1E/file95d242b7da9b/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Thu Feb 19 01:56:10 2026 ─────────────
#> ── Running step: read_realignment @ Thu Feb 19 01:56:10 2026 ───────────────────
#> Realigning sample sample1 -> /tmp/RtmpjRxi1E/file95d242b7da9b/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/RtmpjRxi1E/file95d242b7da9b/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/RtmpjRxi1E/file95d242b7da9b/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Thu Feb 19 01:56:12 2026 ──────────
plot_durations(pipeline)
```
