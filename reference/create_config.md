# Create Configuration File From Arguments

Create Configuration File From Arguments

## Usage

``` r
create_config(outdir, type = "sc_3end", ...)
```

## Arguments

- outdir:

  the destination directory for the configuration file

- type:

  use an example config, available values:

  "sc_3end"

  :   \- config for 10x 3' end ONT reads

  "SIRV"

  :   \- config for the SIRV example reads

- ...:

  Configuration parameters (using dot for nested parameters)

  seed

  :   \- Integer. Seed for minimap2.

  threads

  :   \- Number of threads to use.

  do_barcode_demultiplex

  :   \- Boolean. Specifies whether to run the barcode demultiplexing
      step.

  do_genome_alignment

  :   \- Boolean. Specifies whether to run the genome alignment step.
      `TRUE` is recommended

  do_gene_quantification

  :   \- Boolean. Specifies whether to run gene quantification using the
      genome alignment results. `TRUE` is recommended

  do_isoform_identification

  :   \- Boolean. Specifies whether to run the isoform identification
      step. `TRUE` is recommended

  bambu_isoform_identification

  :   \- Boolean. Whether to use Bambu for isoform identification.

  multithread_isoform_identification

  :   \- Boolean. Whether to use FLAMES' new multithreaded Cpp
      implementation for isoform identification.

  do_read_realignment

  :   \- Boolean. Specifies whether to run the read realignment step.
      `TRUE` is recommended

  do_transcript_quantification

  :   \- Boolean. Specifies whether to run the transcript quantification
      step. `TRUE` is recommended

  barcode_parameters.max_bc_editdistance

  :   \- Maximum edit distance for barcode matching

  barcode_parameters.pattern.primer

  :   \- Primer sequence pattern

  isoform_parameters.max_dist

  :   \- Maximum distance allowed when merging splicing sites

  ...

  :   \- Other nested parameters, using dot to indicate nested section

## Value

file path to the config file created

## Details

Create a list object containing the arguments supplied in a format
usable for the FLAMES pipeline, and writes the object to a JSON file,
which is located with the prefix 'config\_' in the supplied `outdir`.
Default values from `extdata/config_sclr_nanopore_3end.json` will be
used for unprovided parameters.

Parameters can be specified using dot to indicate nested sections, e.g.,
`barcode_parameters.max_bc_editdistance = 3` or
`barcode_parameters.pattern.primer = "ATCG"`. Alternatively, you can
open the created config file and edit it manually.

## Examples

``` r
# create the default configuration file
outdir <- tempdir()
config <- create_config(outdir)
#> Writing configuration parameters to:  /tmp/RtmpaNKtnp/config_file_32976.json 

# create config with custom parameters including nested ones
config <- create_config(outdir,
  threads = 16,
  barcode_parameters.max_bc_editdistance = 3,
  barcode_parameters.pattern.primer = "ATCGATCG",
  isoform_parameters.min_sup_cnt = 10,
  # use the coverage model in oarfish
  # via supplying additional CLI arguments
  additional_arguments.oarfish = c("--model-coverage")
)
#> Writing configuration parameters to:  /tmp/RtmpaNKtnp/config_file_32976.json 
```
