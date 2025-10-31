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
#> Writing configuration parameters to:  /tmp/RtmpM0tt18/filea7097e2e6656/config_file_42761.json 
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- run_FLAMES(pipeline)
#> ── Running step: genome_alignment @ Fri Oct 31 05:41:21 2025 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpM0tt18/filea7097e2e6656/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpM0tt18/filea7097e2e6656/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpM0tt18/filea7097e2e6656/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Fri Oct 31 05:41:21 2025 ─────────────
#> ── Running step: read_realignment @ Fri Oct 31 05:41:21 2025 ───────────────────
#> Realigning sample sample1 -> /tmp/RtmpM0tt18/filea7097e2e6656/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/RtmpM0tt18/filea7097e2e6656/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/RtmpM0tt18/filea7097e2e6656/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Fri Oct 31 05:41:22 2025 ──────────
```
