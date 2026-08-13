# BLAZE Assign reads to cell barcodes.

Uses BLAZE to generate barcode list and assign reads to cell barcodes.

## Usage

``` r
blaze(
  expect_cells,
  fq_in,
  outdir,
  fq_out,
  sample_name = "",
  additional_args = NULL,
  ...
)
```

## Arguments

- expect_cells:

  Integer, expected number of cells. Note: this could be just a rough
  estimate. E.g., the targeted number of cells.

- fq_in:

  File path to the fastq file used as a query sequence file

- outdir:

  Output directory to save BLAZE results.

- fq_out:

  File path to save the output fastq file containing reads assigned to
  cell barcodes.

- sample_name:

  Sample name prefix for output files. Default is an empty string.

- additional_args:

  Additional command line style arguments to be passed to BLAZE. E.g.
  c("–10x-kit-version", "3v3")

- ...:

  Additional BLAZE configuration parameters. E.g., setting
  \`overwrite=TRUE\` is equivalent to switch on the \`–overwrite\`
  option. Note that the specified parameters will override the
  parameters specified in the configuration file. All available options
  can be found at https://github.com/shimlab/BLAZE.

## Value

A `data.frame` summarising the reads aligned. Other outputs are written
to disk. The details of the output files can be found at
https://github.com/shimlab/BLAZE.

## Examples

``` r
outdir <- tempfile()
dir.create(outdir)
fastq <- system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
blaze(
  expect_cells = 10, fastq,
  outdir = outdir,
  fq_out = file.path(outdir, "blaze_matched_reads.fastq.gz"),
  overwrite = TRUE
)
#> $overwrite
#> [1] TRUE
#> 
#> Running BLAZE...
#> Argument:  --expect-cells  10 --output-prefix  /tmp/RtmpgehRKJ/filebbc054663391/ --output-fastq  matched_reads.fastq.gz --overwrite --minimal_stdout   /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz 
#> NULL
```
