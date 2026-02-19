# Get pipeline results

This function returns the results of the pipeline as a
`SummarizedExperiment` object, a `SingleCellExperiment` object, or a
list of `SingleCellExperiment` objects, depending on the pipeline type.

## Usage

``` r
experiment(pipeline)

# S4 method for class 'FLAMES.Pipeline'
experiment(pipeline)
```

## Arguments

- pipeline:

  A FLAMES.Pipeline object.

## Value

A `SummarizedExperiment` object, a `SingleCellExperiment` object, or a
list of `SingleCellExperiment` objects.

## Examples

``` r
pipeline <- example_pipeline(type = "BulkPipeline")
#> Writing configuration parameters to:  /tmp/RtmpjRxi1E/file95d2468277da/config_file_38354.json 
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- run_FLAMES(pipeline)
#> ── Running step: genome_alignment @ Thu Feb 19 01:55:41 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpjRxi1E/file95d2468277da/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpjRxi1E/file95d2468277da/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpjRxi1E/file95d2468277da/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Thu Feb 19 01:55:42 2026 ─────────────
#> ── Running step: read_realignment @ Thu Feb 19 01:55:42 2026 ───────────────────
#> Realigning sample sample1 -> /tmp/RtmpjRxi1E/file95d2468277da/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/RtmpjRxi1E/file95d2468277da/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/RtmpjRxi1E/file95d2468277da/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Thu Feb 19 01:55:44 2026 ──────────
se <- experiment(pipeline)
```
