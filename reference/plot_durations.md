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
#> ℹ Writing configuration to: /tmp/Rtmp4nGYdi/filebc447a45d33f/config_file_48196.json
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- run_FLAMES(pipeline)
#> ── Running step: genome_alignment @ Fri Jun 26 08:24:58 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/Rtmp4nGYdi/filebc447a45d33f/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/Rtmp4nGYdi/filebc447a45d33f/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/Rtmp4nGYdi/filebc447a45d33f/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Fri Jun 26 08:25:00 2026 ─────────────
#> ── Running step: read_realignment @ Fri Jun 26 08:25:00 2026 ───────────────────
#> Realigning sample sample1 -> /tmp/Rtmp4nGYdi/filebc447a45d33f/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/Rtmp4nGYdi/filebc447a45d33f/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/Rtmp4nGYdi/filebc447a45d33f/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Fri Jun 26 08:25:01 2026 ──────────
plot_durations(pipeline)
```
