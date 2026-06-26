# Filter transcript coverage

Filter the transcript coverage by applying a filter function to the
coverage values.

## Usage

``` r
filter_coverage(x, filter_fn = convolution_filter)
```

## Arguments

- x:

  The tibble returned by
  [`get_coverage`](https://mritchielab.github.io/FLAMES/reference/get_coverage.md),
  or a BAM file path, or a GAlignments object.

- filter_fn:

  The filter function to apply to the coverage values. The function
  should take a numeric vector of coverage values and return a logical
  value (TRUE if the transcript passes the filter, FALSE otherwise). The
  default filter function is
  [`convolution_filter`](https://mritchielab.github.io/FLAMES/reference/convolution_filter.md),
  which filters out transcripts with sharp drops / rises in coverage.

## Value

a tibble of the transcript information and coverages, with transcipts
that pass the filter

## Examples

``` r
ppl <- example_pipeline("BulkPipeline")
#> ℹ Writing configuration to: /tmp/Rtmp4nGYdi/filebc4424583a0f/config_file_48196.json
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
steps(ppl)["isoform_identification"] <- FALSE
ppl <- run_step(ppl, "read_realignment")
#> ── Running step: read_realignment @ Fri Jun 26 08:24:30 2026 ───────────────────
#> Using reference annotation for transcriptome assembly.
#> Realigning sample sample1 -> /tmp/Rtmp4nGYdi/filebc4424583a0f/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/Rtmp4nGYdi/filebc4424583a0f/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/Rtmp4nGYdi/filebc4424583a0f/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
x <- get_coverage(ppl@transcriptome_bam[[1]])
nrow(x)
#> [1] 2
filter_coverage(x) |>
  nrow()
#> 2 transcripts found in the BAM file.
#> 0(0%) transcripts failed the filter.
#> Failed transcripts account for 0 reads, out of 219(0%) reads in total.
#> [1] 2
```
