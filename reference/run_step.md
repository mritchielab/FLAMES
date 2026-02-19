# Execute a single step of the FLAMES pipeline

This function runs the specified step of the FLAMES pipeline.

## Usage

``` r
run_step(pipeline, step, disable_controller = TRUE)

# S4 method for class 'FLAMES.Pipeline'
run_step(pipeline, step, disable_controller = TRUE)
```

## Arguments

- pipeline:

  A FLAMES.Pipeline object.

- step:

  The step to run. One of "barcode_demultiplex", "genome_alignment",
  "gene_quantification", "isoform_identification", "read_realignment",
  or "transcript_quantification".

- disable_controller:

  (optional) If TRUE, the step will be executed in the current R
  session, instead of using crew controllers.

## Value

An updated FLAMES.Pipeline object.

## See also

[`run_FLAMES`](https://mritchielab.github.io/FLAMES/reference/run_FLAMES.md)
to run the entire pipeline.
[`resume_FLAMES`](https://mritchielab.github.io/FLAMES/reference/resume_FLAMES.md)
to resume a pipeline from the last completed step.

## Examples

``` r
pipeline <- example_pipeline("BulkPipeline")
#> Writing configuration parameters to:  /tmp/RtmpjRxi1E/file95d23e4ba071/config_file_38354.json 
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- run_step(pipeline, "genome_alignment")
#> ── Running step: genome_alignment @ Thu Feb 19 01:56:44 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpjRxi1E/file95d23e4ba071/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpjRxi1E/file95d23e4ba071/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpjRxi1E/file95d23e4ba071/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
```
