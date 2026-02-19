# Execute a FLAMES pipeline

This function runs the FLAMES pipeline. It will run all steps in the
pipeline.

## Usage

``` r
run_FLAMES(pipeline, overwrite = FALSE)

# S4 method for class 'FLAMES.Pipeline'
run_FLAMES(pipeline, overwrite = FALSE)
```

## Arguments

- pipeline:

  A FLAMES.Pipeline object.

- overwrite:

  (optional) If TRUE, the pipeline will be re-run even if some steps are
  already completed.

## Value

An updated FLAMES.Pipeline object.

## See also

[`resume_FLAMES`](https://mritchielab.github.io/FLAMES/reference/resume_FLAMES.md)
to resume a pipeline from the last completed step.

## Examples

``` r
pipeline <- example_pipeline("BulkPipeline")
#> Writing configuration parameters to:  /tmp/RtmpjRxi1E/file95d278147a4d/config_file_38354.json 
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- run_FLAMES(pipeline)
#> ── Running step: genome_alignment @ Thu Feb 19 01:56:41 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpjRxi1E/file95d278147a4d/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpjRxi1E/file95d278147a4d/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpjRxi1E/file95d278147a4d/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Thu Feb 19 01:56:42 2026 ─────────────
#> ── Running step: read_realignment @ Thu Feb 19 01:56:42 2026 ───────────────────
#> Realigning sample sample1 -> /tmp/RtmpjRxi1E/file95d278147a4d/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/RtmpjRxi1E/file95d278147a4d/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/RtmpjRxi1E/file95d278147a4d/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Thu Feb 19 01:56:44 2026 ──────────
```
