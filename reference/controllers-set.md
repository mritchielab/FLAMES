# Set controllers

Sets the controllers for the pipeline.

## Usage

``` r
controllers(pipeline) <- value

# S4 method for class 'FLAMES.Pipeline'
controllers(pipeline) <- value
```

## Arguments

- pipeline:

  A FLAMES.Pipeline object.

- value:

  A `crew_class_controller` object or a named list of
  `crew_class_controller` objects. If a single controller is provided,
  it will be used for all steps in the pipeline. If a named list is
  provided, steps with names that match the names of the list will use
  the corresponding controller, and steps without a specified controller
  will use the current R session.

## Value

An updated FLAMES.Pipeline object with the specified controllers.

## Examples

``` r
pipeline <- example_pipeline()
#> ℹ Writing configuration to: /tmp/RtmpgehRKJ/filebbc037a4f742/config_file_48064.json
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
# Only set the genome alignment controller
controllers(pipeline) <- list(genome_alignment = crew::crew_controller_local())
# Same as above
controllers(pipeline)[["genome_alignment"]] <- crew::crew_controller_local()
# Set a controller for all steps
controllers(pipeline) <- crew::crew_controller_local()
# Unset all controllers and use the current R session
controllers(pipeline) <- list()
```
