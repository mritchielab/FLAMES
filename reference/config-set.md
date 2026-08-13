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
#> ℹ Writing configuration to: /tmp/RtmpgehRKJ/filebbc04b1a8824/config_file_48064.json
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
# Set a new configuration
config(pipeline) <- create_config(outdir = tempdir())
#> ℹ Writing configuration to: /tmp/RtmpgehRKJ/config_file_48064.json
```
