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
#> Writing configuration parameters to:  /tmp/RtmpM0tt18/filea709637feacf/config_file_42761.json 
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- run_FLAMES(pipeline)
#> ── Running step: genome_alignment @ Fri Oct 31 05:40:50 2025 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpM0tt18/filea709637feacf/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpM0tt18/filea709637feacf/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpM0tt18/filea709637feacf/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Fri Oct 31 05:40:51 2025 ─────────────
#> ── Running step: read_realignment @ Fri Oct 31 05:40:51 2025 ───────────────────
#> Realigning sample sample1 -> /tmp/RtmpM0tt18/filea709637feacf/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/RtmpM0tt18/filea709637feacf/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/RtmpM0tt18/filea709637feacf/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Fri Oct 31 05:40:51 2025 ──────────
plot_durations(pipeline)
```
