# Transcript quantification

Calculate the transcript count matrix by parsing the re-alignment file.

## Usage

``` r
quantify_transcript(
  annotation,
  outdir,
  config,
  pipeline = "sc_single_sample",
  ...
)
```

## Arguments

- annotation:

  The file path to the annotation file in GFF3 format

- outdir:

  The path to directory to store all output files.

- config:

  Parsed FLAMES configurations.

- pipeline:

  The pipeline type as a character string, either `sc_single_sample`
  (single-cell, single-sample),

- ...:

  Supply sample names as character vector (e.g.
  `samples = c("name1", "name2", ...)`) for muti-sample or bulk
  pipeline. `bulk` (bulk, single or multi-sample), or `sc_multi_sample`
  (single-cell, multiple samples)

## Value

A `SingleCellExperiment` object for single-cell pipeline, a list of
`SingleCellExperiment` objects for multi-sample pipeline, or a
`SummarizedExperiment` object for bulk pipeline.
