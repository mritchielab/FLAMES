# Index the reference genome for minimap2

Calls minimap2 to index the reference genome.

## Usage

``` r
index_genome(pipeline, path, additional_args = c("-k", "14"))

# S4 method for class 'FLAMES.Pipeline'
index_genome(pipeline, path, additional_args = c("-k", "14"))
```

## Arguments

- pipeline:

  A FLAMES.Pipeline object.

- path:

  The file path to save the minimap2 index. If not provided, it will be
  saved to the output directory with the name "genome.mmi".

- additional_args:

  (optional) Additional arguments to pass to minimap2.

## Value

A `SummarizedExperiment` object, a `SingleCellExperiment` object, or a
list of `SingleCellExperiment` objects.

## Examples

``` r
pipeline <- example_pipeline(type = "BulkPipeline")
#> Writing configuration parameters to:  /tmp/RtmpaNKtnp/file80d06fef1b83/config_file_32976.json 
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline <- index_genome(pipeline)
```
