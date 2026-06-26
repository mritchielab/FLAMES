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
#> ℹ Writing configuration to: /tmp/Rtmp4nGYdi/filebc4435a12ce/config_file_48196.json
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- run_FLAMES(pipeline)
#> ── Running step: genome_alignment @ Fri Jun 26 08:25:27 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/Rtmp4nGYdi/filebc4435a12ce/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/Rtmp4nGYdi/filebc4435a12ce/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/Rtmp4nGYdi/filebc4435a12ce/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Fri Jun 26 08:25:28 2026 ─────────────
#> ── Running step: read_realignment @ Fri Jun 26 08:25:28 2026 ───────────────────
#> Realigning sample sample1 -> /tmp/Rtmp4nGYdi/filebc4435a12ce/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/Rtmp4nGYdi/filebc4435a12ce/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/Rtmp4nGYdi/filebc4435a12ce/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Fri Jun 26 08:25:30 2026 ──────────
```
