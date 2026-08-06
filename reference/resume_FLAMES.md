# Resume a FLAMES pipeline

This function resumes a FLAMES pipeline by running configured but
unfinished steps.

## Usage

``` r
resume_FLAMES(pipeline)

# S4 method for class 'FLAMES.Pipeline'
resume_FLAMES(pipeline)
```

## Arguments

- pipeline:

  A FLAMES.Pipeline object.

## Value

An updated FLAMES.Pipeline object.

## See also

[`run_FLAMES`](https://mritchielab.github.io/FLAMES/reference/run_FLAMES.md)
to run the entire pipeline.

## Examples

``` r
pipeline <- example_pipeline("BulkPipeline")
#> ℹ Writing configuration to: /tmp/RtmpnC89xy/filebc155fe29f59/config_file_48149.json
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- run_step(pipeline, "genome_alignment")
#> ── Running step: genome_alignment @ Thu Aug  6 04:49:34 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpnC89xy/filebc155fe29f59/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpnC89xy/filebc155fe29f59/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpnC89xy/filebc155fe29f59/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
pipeline <- resume_FLAMES(pipeline)
#> Resuming pipeline from step: isoform_identification
#> ── Running step: isoform_identification @ Thu Aug  6 04:49:35 2026 ─────────────
#> ── Running step: read_realignment @ Thu Aug  6 04:49:35 2026 ───────────────────
#> Realigning sample sample1 -> /tmp/RtmpnC89xy/filebc155fe29f59/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/RtmpnC89xy/filebc155fe29f59/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/RtmpnC89xy/filebc155fe29f59/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Thu Aug  6 04:49:37 2026 ──────────
```
