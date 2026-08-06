# Calculate mutation positions within the gene body

Given a set of mutations and gene annotation, calculate the position of
each mutation within the gene body they are in.

## Usage

``` r
mutation_positions(
  mutations,
  annotation,
  type = "relative",
  bin = FALSE,
  by = c(region = "gene_name"),
  threads = 1
)
```

## Arguments

- mutations:

  either the tibble output from `find_variants`. It must have columns
  `seqnames`, `pos`, and a third column for specifying the gene id or
  gene name. The mutation must be within the gene region.

- annotation:

  Either path to the annotation file (GTF/GFF) or a GRanges object of
  the gene annotation.

- type:

  character(1): the type of position to calculate. Can be one of "TSS"
  (distance from the transcription start site), "TES" (distance from the
  transcription end site), or "relative" (relative position within the
  gene body).

- bin:

  logical(1): whether to bin the relative positions into 100 bins. Only
  applicable when `type = "relative"`.

- by:

  character(1): the column name in the annotation to match with the gene
  annotation. E.g. `c("region" = "gene_name")` to match the \`region\`
  column in the mutations with the \`gene_name\` column in the
  annotation.

- threads:

  integer(1): number of threads to use.

## Value

A numeric vector of positions of each mutation within the gene body.
When `type = "relative"`, the positions are normalized to the gene
length, ranging from 0 (start of the gene) to 1 (end of the gene). When
`type = "TSS"` / `type = "TES"`, the distances from the transcription
start / end site. If `bin = TRUE`, and `type = "relative"`, the relative
positions are binned into 100 bins along the gene body, and the output
is a matrix with the number of mutations in each bin, the rows are named
by the `by` column (e.g. gene name).

## Examples

``` r
variants <- data.frame(
  seqnames = rep("chr14", 8),
  pos = c(1084, 1085, 1217, 1384, 2724, 2789, 5083, 5147),
  region = rep("Rps24", 8)
)
positions <-
 mutation_positions(
   mutations = variants,
   annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES")
 )
#> 04:49:01 Reading annotation ...
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#> 
```
