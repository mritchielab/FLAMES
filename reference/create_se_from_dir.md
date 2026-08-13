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
#> ℹ Writing configuration to: /tmp/RtmpgehRKJ/filebbc07978e040/config_file_48064.json
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
ppl <- run_FLAMES(ppl)
#> FLAMES version 2.7.1 (unknown source)
#> ── Running step: genome_alignment @ Thu Aug 13 10:09:39 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample sample1 -> /tmp/RtmpgehRKJ/filebbc07978e040/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample2 -> /tmp/RtmpgehRKJ/filebbc07978e040/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample sample3 -> /tmp/RtmpgehRKJ/filebbc07978e040/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Thu Aug 13 10:09:40 2026 ─────────────
#> ── Running step: read_realignment @ Thu Aug 13 10:09:40 2026 ───────────────────
#> Realigning sample sample1 -> /tmp/RtmpgehRKJ/filebbc07978e040/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/RtmpgehRKJ/filebbc07978e040/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/RtmpgehRKJ/filebbc07978e040/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> ── Running step: transcript_quantification @ Thu Aug 13 10:09:42 2026 ──────────
se1 <- experiment(ppl)
se2 <- create_se_from_dir(ppl@outdir, ppl@annotation)
```
