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
#> Writing configuration parameters to:  /tmp/RtmpmJ8vO7/file80d052bdbc9a/config_file_32976.json 
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- run_step(pipeline, "genome_alignment")
#> ── Running step: genome_alignment @ Fri Oct 31 06:52:23 2025 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpmJ8vO7/file80d052bdbc9a/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpmJ8vO7/file80d052bdbc9a/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpmJ8vO7/file80d052bdbc9a/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
pipeline <- resume_FLAMES(pipeline)
#> Resuming pipeline from step: isoform_identification
#> ── Running step: isoform_identification @ Fri Oct 31 06:52:24 2025 ─────────────
#> ── Running step: read_realignment @ Fri Oct 31 06:52:24 2025 ───────────────────
#> Realigning sample sample1 -> /tmp/RtmpmJ8vO7/file80d052bdbc9a/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/RtmpmJ8vO7/file80d052bdbc9a/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/RtmpmJ8vO7/file80d052bdbc9a/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Fri Oct 31 06:52:24 2025 ──────────
```
