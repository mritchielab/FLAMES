# BLAZE Assign reads to cell barcodes.

Uses BLAZE to generate barcode list and assign reads to cell barcodes.

## Usage

``` r
blaze(expect_cells, fq_in, additional_args = NULL, ...)
```

## Arguments

- expect_cells:

  Integer, expected number of cells. Note: this could be just a rough
  estimate. E.g., the targeted number of cells.

- fq_in:

  File path to the fastq file used as a query sequence file

- additional_args:

  Additional command line style arguments to be passed to BLAZE. E.g.
  c("–10x-kit-version", "3v3")

- ...:

  Additional BLAZE configuration parameters. E.g., setting
  \`'output-prefix'='some_prefix'\` is equivalent to specifying
  \`–output-prefix some_prefix\` in BLAZE; Similarly, \`overwrite=TRUE\`
  is equivalent to switch on the \`–overwrite\` option. Note that the
  specified parameters will override the parameters specified in the
  configuration file. All available options can be found at
  https://github.com/shimlab/BLAZE.

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
  "output-prefix" = file.path(outdir, ""),
  "output-fastq" = file.path(outdir, "output.fastq"),
  overwrite=TRUE
)
#> $`output-prefix`
#> [1] "/tmp/RtmpM0tt18/filea7097b4b7c56/"
#> 
#> $`output-fastq`
#> [1] "/tmp/RtmpM0tt18/filea7097b4b7c56/output.fastq"
#> 
#> $overwrite
#> [1] TRUE
#> 
#> Running BLAZE...
#> Argument:  --expect-cells  10 --overwrite --minimal_stdout  --output-prefix /tmp/RtmpM0tt18/filea7097b4b7c56/ --output-fastq /tmp/RtmpM0tt18/filea7097b4b7c56/output.fastq  /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz 
```
