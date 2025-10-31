# mutation positions within the gene body

Given a set of mutations and a gene annotation, calculate the position
of each mutation within the gene body. The gene annotation must have the
following types: "gene" and "exon". The gene annotation must be for one
gene only. The mutations must be within the gene region. The function
will merge overlapping exons and calculate the position of each mutation
within the gene body, excluding intronic regions.

## Usage

``` r
mutation_positions_single(mutations, annotation_grange, type, verbose = TRUE)
```

## Arguments

- mutations:

  either the tibble output from `find_variants` or a GRanges object.
  Make sure to filter it for only the gene of interest.

- annotation_grange:

  GRanges: the gene annotation. Must have the following types: "gene"
  and "exon".

- type:

  character(1): the type of position to calculate. Can be one of "TSS"
  (distance from the transcription start site), "TES" (distance from the
  transcription end site), or "relative" (relative position within the
  gene body).

- verbose:

  logical(1): whether to print messages.

## Value

A numeric vector of positions of each mutation within the gene body.
When `type = "relative"`, the positions are normalized to the gene
length, ranging from 0 (start of the gene) to 1 (end of the gene). When
`type = "TSS"` / `type = "TES"`, the distances from the transcription
start / end site.
