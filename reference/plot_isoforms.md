# Plot isoforms

Plot isoforms, either from a gene or a list of transcript ids.

## Usage

``` r
plot_isoforms(
  sce,
  gene_id,
  transcript_ids,
  n = 4,
  format = "plot_grid",
  colors
)
```

## Arguments

- sce:

  The `SingleCellExperiment` object containing transcript counts,
  `rowRanges` and `rowData` with `gene_id` and `transcript_id` columns.

- gene_id:

  The gene symbol of interest, ignored if `transcript_ids` is provided.

- transcript_ids:

  The transcript ids to plot.

- n:

  The number of top isoforms to plot from the gene. Ignored if
  `transcript_ids` is provided.

- format:

  The format of the output, either "plot_grid" or "list".

- colors:

  A character vector of colors to use for the isoforms. If not provided,
  gray will be used. for all isoforms.

## Value

When `format = "list"`, a list of `ggplot` objects is returned.
Otherwise, a grid of the plots is returned.

## Details

This function takes a `SingleCellExperiment` object and plots the top
isoforms of a gene, or a list of specified transcript ids. Either as a
list of plots or together in a grid. This function wraps the
[`ggbio::geom_alignment`](https://rdrr.io/pkg/ggbio/man/geom_alignment-method.html)
function to plot the isoforms, and orders the isoforms by expression
levels (when specifying a gene) or by the order of the transcript_ids.

## Examples

``` r
data(scmixology_lib10_transcripts)
plot_isoforms(scmixology_lib10_transcripts, gene_id = "ENSG00000108107")
#> Constructing graphics...
#> Constructing graphics...
#> Constructing graphics...
#> Constructing graphics...

```
