# Get read coverages from BAM file

Get the read coverages for each transcript in the BAM file (or a
GAlignments object). The read coverages are sampled at 100 positions
along the transcript, and the coverage is scaled by dividing the
coverage at each position by the total read counts for the transcript.
If a BAM file is provided, alignment with MAPQ \< 5, secondary
alignments and supplementary alignments are filtered out. A
`GAlignments` object can also be provided in case alternative filtering
is desired.

## Usage

``` r
get_coverage(bam, min_counts = 10, remove_UTR = FALSE, annotation)
```

## Arguments

- bam:

  path to the BAM file, or a parsed GAlignments object

- min_counts:

  numeric, the minimum number of alignments required for a transcript to
  be included

- remove_UTR:

  logical, if `TRUE`, remove the UTRs from the coverage

- annotation:

  (Required if `remove_UTR = TRUE`) path to the GTF annotation file

## Value

a tibble of the transcript information and coverages, with the following
columns:

- `transcript`: the transcript name / ID

- `read_counts`: the total number of aligments for the transcript

- `coverage_1-100`: the coverage at each of the 100 positions along the
  transcript

- `tr_length`: the length of the transcript

## Examples

``` r
ppl <- example_pipeline("BulkPipeline")
#> Writing configuration parameters to:  /tmp/Rtmpbzssfl/filea709223994/config_file_42761.json 
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
steps(ppl)["isoform_identification"] <- FALSE
ppl <- run_step(ppl, "read_realignment")
#> ── Running step: read_realignment @ Fri Oct 31 06:02:04 2025 ───────────────────
#> Using reference annotation for transcriptome assembly.
#> Realigning sample sample1 -> /tmp/Rtmpbzssfl/filea709223994/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample2 -> /tmp/Rtmpbzssfl/filea709223994/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
#> Realigning sample sample3 -> /tmp/Rtmpbzssfl/filea709223994/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Skipped sorting BAM files.
x <- get_coverage(ppl@transcriptome_bam[[1]])
head(x)
#> # A tibble: 1 × 103
#>   transcript   coverage_1 coverage_2 coverage_3 coverage_4 coverage_5 coverage_6
#>   <chr>             <dbl>      <dbl>      <dbl>      <dbl>      <dbl>      <dbl>
#> 1 ENSMUST0000…      0.880          1          1          1          1          1
#> # ℹ 96 more variables: coverage_7 <dbl>, coverage_8 <dbl>, coverage_9 <dbl>,
#> #   coverage_10 <dbl>, coverage_11 <dbl>, coverage_12 <dbl>, coverage_13 <dbl>,
#> #   coverage_14 <dbl>, coverage_15 <dbl>, coverage_16 <dbl>, coverage_17 <dbl>,
#> #   coverage_18 <dbl>, coverage_19 <dbl>, coverage_20 <dbl>, coverage_21 <dbl>,
#> #   coverage_22 <dbl>, coverage_23 <dbl>, coverage_24 <dbl>, coverage_25 <dbl>,
#> #   coverage_26 <dbl>, coverage_27 <dbl>, coverage_28 <dbl>, coverage_29 <dbl>,
#> #   coverage_30 <dbl>, coverage_31 <dbl>, coverage_32 <dbl>, …
```
