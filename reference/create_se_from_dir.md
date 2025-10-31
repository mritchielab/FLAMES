# Create `SummarizedExperiment` object from `FLAMES` output folder

Create `SummarizedExperiment` object from `FLAMES` output folder

## Usage

``` r
create_se_from_dir(outdir, annotation, quantification = "FLAMES")
```

## Arguments

- outdir:

  The folder containing `FLAMES` output files

- annotation:

  (Optional) the annotation file that was used to produce the output
  files

- quantification:

  (Optional) the quantification method used to generate the output files
  (either "FLAMES" or "Oarfish".). If not specified, the function will
  attempt to determine the quantification method.

## Value

a `SummarizedExperiment` object

## Examples

``` r
ppl <- example_pipeline("BulkPipeline")
#> Writing configuration parameters to:  /tmp/RtmpM0tt18/filea7093c0ab44/config_file_42761.json 
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
ppl <- run_FLAMES(ppl)
#> ── Running step: genome_alignment @ Fri Oct 31 05:40:20 2025 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpM0tt18/filea7093c0ab44/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpM0tt18/filea7093c0ab44/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpM0tt18/filea7093c0ab44/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Fri Oct 31 05:40:20 2025 ─────────────
#> ── Running step: read_realignment @ Fri Oct 31 05:40:21 2025 ───────────────────
#> Realigning sample sample1 -> /tmp/RtmpM0tt18/filea7093c0ab44/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/RtmpM0tt18/filea7093c0ab44/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/RtmpM0tt18/filea7093c0ab44/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Fri Oct 31 05:40:21 2025 ──────────
se1 <- experiment(ppl)
se2 <- create_se_from_dir(ppl@outdir, ppl@annotation)
```
