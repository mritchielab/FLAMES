# FLAMES isoform reduced dimensions plots

Plot expression of top n isoforms of a gene in reduced dimensions

## Usage

``` r
plot_isoform_reduced_dim(
  sce,
  gene_id,
  transcript_ids,
  n = 4,
  reduced_dim_name = "UMAP",
  use_gene_dimred = FALSE,
  expr_func = function(x) {
     SingleCellExperiment::logcounts(x)
 },
  col_low = "#313695",
  col_mid = "#FFFFBF",
  col_high = "#A50026",
  alpha = 0.5,
  size = 0.2,
  ggtheme = theme_minimal() + theme(axis.text = element_blank()),
  color_quantile = 1,
  format = "plot_grid",
  ...
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

- reduced_dim_name:

  The name of the reduced dimension to use for plotting cells.

- use_gene_dimred:

  Whether to use gene-level reduced dimensions for plotting. Set to
  `TRUE` if the `SingleCellExperiment` has gene counts in main assay and
  transcript counts in `altExp`.

- expr_func:

  The function to extract expression values from the
  `SingleCellExperiment` object. Default is `logcounts`. Alternatively,
  `counts` can be used for raw counts.

- col_low:

  Color for cells with low expression levels in UMAPs.

- col_mid:

  Color for cells with intermediate expression levels in UMAPs.

- col_high:

  Color for cells with high expression levels in UMAPs.

- alpha:

  The transparency of the points in the UMAPs.

- size:

  The size of the points in the UMAPs.

- ggtheme:

  The theme to use for the UMAPs.

- color_quantile:

  The lower and upper expression quantile to be displayed bewteen
  `col_low` and `col_high`, e.g. with `color_quantile = 0.95`, cells
  with expressions higher than 95% of other cells will all be shown in
  `col_high`, and cells with expression lower than 95% of other cells
  will all be shown in `col_low`.

- format:

  The format of the output, either "plot_grid" or "list".

- ...:

  Additional arguments to pass to `plot_grid`.

## Value

a `ggplot` object of the UMAP(s)

## Details

Takes `SingleCellExperiment` object and plots an expression on reduced
dimensions with the isoform visualizations along genomic coordinates.

## Examples

``` r
data(scmixology_lib10_transcripts, scmixology_lib10, scmixology_lib90)
scmixology_lib10 <- 
  scmixology_lib10[, colSums(SingleCellExperiment::counts(scmixology_lib10)) > 0]
sce_lr <- scmixology_lib10[, colnames(scmixology_lib10) %in% colnames(scmixology_lib10_transcripts)]
SingleCellExperiment::altExp(sce_lr, "transcript") <-
  scmixology_lib10_transcripts[, colnames(sce_lr)]
combined_sce <- combine_sce(sce_lr, scmixology_lib90)
combined_sce <- combined_sce |>
  scrapper::normalizeRnaCounts.se() |>
  scater::runPCA() |>
  scater::runUMAP()
#> using unknown matrix fallback for 'data.frame'
#> using unknown matrix fallback for 'data.frame'
combined_imputed_sce <- sc_impute_transcript(combined_sce)
#> Warning: 'buildSNNGraph' is deprecated.
#> Use 'bluster::makeSNNGraph' instead.
#> See help("Deprecated")
#> Imputing transcript counts ...
plot_isoform_reduced_dim(combined_sce, 'ENSG00000108107')

plot_isoform_reduced_dim(combined_imputed_sce, 'ENSG00000108107')

```
