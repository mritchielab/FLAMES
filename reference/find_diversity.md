# Compute Gene Isoform Entropy Matrix

Calculates normalized Shannon entropy for gene isoform expression across
cells. Higher entropy indicates more diverse isoform usage, lower
entropy indicates dominance by fewer isoforms.

## Usage

``` r
find_diversity(
  sce,
  assay = "counts",
  gene_col = "gene_id",
  alpha = .Machine$double.xmin,
  min_counts_per_cell = 5,
  isoform_min_pct_cells = 0.05,
  isoform_cumulative_pct = 0.95,
  min_cell_fraction = 0.25,
  threads = 1,
  show_progress = interactive()
)
```

## Arguments

- sce:

  A `SingleCellExperiment` object

- assay:

  Name of assay containing isoform counts (default: "counts")

- gene_col:

  Column name in rowData containing gene identifiers (default:
  "gene_id")

- alpha:

  Pseudocount added to avoid log(0) (default: .Machine\$double.xmin)

- min_counts_per_cell:

  Minimum total gene counts per cell to include (default: 5)

- isoform_min_pct_cells:

  Minimum fraction of cells expressing each isoform (default: 0.05)

- isoform_cumulative_pct:

  Keep top isoforms contributing to this cumulative proportion (default:
  0.95)

- min_cell_fraction:

  Minimum fraction of cells with valid entropy per gene (default: 0.25)

- threads:

  Number of threads for parallel processing (default: 1)

- show_progress:

  Logical indicating whether to show progress (default: TRUE if
  interactive)

## Value

Matrix with genes as rows and cells as columns containing normalized
entropy values (0-1).

## Examples

``` r
sce <- scuttle::mockSCE(ncells = 50, ngenes = 30)
SummarizedExperiment::rowData(sce)$gene_id <- sort(
  paste0("gene", sample(1:9, nrow(sce), replace = TRUE))
)
res <- find_diversity(sce, threads = 2)
```
