# Add gene counts to a `SingleCellExperiment` object

Add gene counts to a `SingleCellExperiment` object as an `altExps` slot
named `gene`.

## Usage

``` r
add_gene_counts(sce, gene_count_file)
```

## Arguments

- sce:

  A `SingleCellExperiment` object.

- gene_count_file:

  The file path to the gene count file. If missing, the function will
  try to find the gene count file in the output directory.

## Value

A `SingleCellExperiment` object with gene counts added.

## Examples

``` r
# Set up a mock SingleCellExperiment object
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = matrix(0, nrow = 10, ncol = 10))
)
colnames(sce) <- paste0("cell", 1:10)
# Set up a mock gene count file
gene_count_file <- tempfile()
gene_mtx <- matrix(1:10, nrow = 2, ncol = 5)
colnames(gene_mtx) <- paste0("cell", 1:5)
rownames(gene_mtx) <- c("gene1", "gene2")
write.csv(gene_mtx, gene_count_file)
# Add gene counts to the SingleCellExperiment object
sce <- add_gene_counts(sce, gene_count_file)
# verify the gene counts are added
SingleCellExperiment::altExps(sce)$gene
#> class: SingleCellExperiment 
#> dim: 2 10 
#> metadata(0):
#> assays(1): counts
#> rownames(2): gene1 gene2
#> rowData names(0):
#> colnames(10): cell1 cell2 ... cell9 cell10
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```
