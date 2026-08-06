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
#> ℹ Writing configuration to: /tmp/RtmpnC89xy/filebc1571d490e/config_file_48149.json
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
steps(ppl)
#>       barcode_demultiplex          genome_alignment       gene_quantification 
#>                      TRUE                      TRUE                      TRUE 
#>    isoform_identification          read_realignment transcript_quantification 
#>                      TRUE                      TRUE                      TRUE 
```
