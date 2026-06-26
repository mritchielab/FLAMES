# Add gene counts to a `SingleCellExperiment` object

Add gene counts to a `SingleCellExperiment` object as an `altExps` slot
named `gene`.

## Usage

``` r
add_gene_counts(sce, gene_count_stem)
```

## Arguments

- sce:

  A `SingleCellExperiment` object.

- gene_count_stem:

  The file path stem for the gene count MTX files. The function expects
  three files: `<gene_count_stem>.mtx` (count matrix in Matrix Market
  format), `<gene_count_stem>_features.tsv` (gene names, one per line),
  and `<gene_count_stem>_barcodes.tsv` (barcode names, one per line).

## Value

A `SingleCellExperiment` object with gene counts added as
`altExps(sce)$gene`.

## Examples

``` r
# Set up a mock SingleCellExperiment object
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = matrix(0, nrow = 10, ncol = 10))
)
colnames(sce) <- paste0("cell", 1:10)
# Write mock gene count MTX files
gene_count_stem <- file.path(tempdir(), "gene_count")
gene_mtx <- Matrix::Matrix(1:10, nrow = 2, ncol = 5, sparse = TRUE)
colnames(gene_mtx) <- paste0("cell", 1:5)
rownames(gene_mtx) <- c("gene1", "gene2")
Matrix::writeMM(gene_mtx, paste0(gene_count_stem, ".mtx"))
#> NULL
writeLines(rownames(gene_mtx), paste0(gene_count_stem, "_features.tsv"))
writeLines(colnames(gene_mtx), paste0(gene_count_stem, "_barcodes.tsv"))
# Add gene counts to the SingleCellExperiment object
sce <- add_gene_counts(sce, gene_count_stem)
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
