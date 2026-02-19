# Steps to perform in the pipeline

Steps to perform in the pipeline

## Usage

``` r
steps(pipeline)

# S4 method for class 'FLAMES.Pipeline'
steps(pipeline)
```

## Arguments

- pipeline:

  An object of class \`FLAMES.Pipeline\`

## Value

A named logical vector containing all possible steps for the pipeline.
The names of the vector are the step names, and the values are logical
indicating whether the step is configured to be performed.

## Examples

``` r
ppl <- example_pipeline()
#> Writing configuration parameters to:  /tmp/RtmpjRxi1E/file95d234cffe78/config_file_38354.json 
#> Warning: You have set to use oarfish quantification without gene quantification. Oarfish currently does not collapse UMIs, and gene quantification performs UMI collapsing. You may want to set do_gene_quantification to TRUE for more accurate results.
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: FALSE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
steps(ppl)
#>       barcode_demultiplex          genome_alignment       gene_quantification 
#>                      TRUE                      TRUE                     FALSE 
#>    isoform_identification          read_realignment transcript_quantification 
#>                      TRUE                      TRUE                      TRUE 
```
