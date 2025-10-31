# Set pipeline configurations

This function sets the configuration of the pipeline.

## Usage

``` r
config(pipeline) <- value

# S4 method for class 'FLAMES.Pipeline'
config(pipeline) <- value
```

## Arguments

- pipeline:

  An pipeline of class \`FLAMES.Pipeline\`.

- value:

  A list containing the configuration of the pipeline, or a path to a
  JSON configuration file.

## Value

An pipeline of class \`FLAMES.Pipeline\` with the updated configuration.

## Examples

``` r
pipeline <- example_pipeline(type = "BulkPipeline")
#> Writing configuration parameters to:  /tmp/RtmpaNKtnp/file80d012a94631/config_file_32976.json 
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
# Set a new configuration
config(pipeline) <- create_config(outdir = tempdir())
#> Writing configuration parameters to:  /tmp/RtmpaNKtnp/config_file_32976.json 
```
