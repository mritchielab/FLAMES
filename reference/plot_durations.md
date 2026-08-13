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
#> ℹ Writing configuration to: /tmp/RtmpgehRKJ/filebbc02f539113/config_file_48064.json
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- run_FLAMES(pipeline)
#> FLAMES version 2.7.1 (unknown source)
#> ── Running step: genome_alignment @ Thu Aug 13 10:10:17 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpgehRKJ/filebbc02f539113/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpgehRKJ/filebbc02f539113/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpgehRKJ/filebbc02f539113/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Thu Aug 13 10:10:18 2026 ─────────────
#> ── Running step: read_realignment @ Thu Aug 13 10:10:18 2026 ───────────────────
#> Realigning sample sample1 -> /tmp/RtmpgehRKJ/filebbc02f539113/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/RtmpgehRKJ/filebbc02f539113/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/RtmpgehRKJ/filebbc02f539113/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Thu Aug 13 10:10:20 2026 ──────────
plot_durations(pipeline)
```
