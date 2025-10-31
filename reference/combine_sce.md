# Combine SCE

Combine FLT-seq SingleCellExperiment objects

## Usage

``` r
combine_sce(sce_with_lr, sce_without_lr)
```

## Arguments

- sce_with_lr:

  A `SingleCellExperiment` object with both long and short reads. The
  long-read transcript counts should be stored in the 'transcript'
  altExp slot.

- sce_without_lr:

  A `SingleCellExperiment` object with only short reads.

## Value

A `SingleCellExperiment` object with combined gene counts and a
"transcript" altExp slot.

## Details

For protcols like FLT-seq that generate two libraries, one with both
short and long reads, and one with only short reads, this function
combines the two libraries into a single `SingleCellExperiment` object.
For the library with both long and short reads, the long-read transcript
counts should be stored in the 'transcript' altExp slot of the
`SingleCellExperiment` object. This function will combine the short-read
gene counts of both libraries, and for the transcripts counts, it will
leave `NA` values for the cells from the short-read only library. The
`sc_impute_transcript` function can then be used to impute the `NA`
values.

## Examples

``` r
with_lr <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(rpois(100, 5), ncol = 10)))
without_lr <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(rpois(200, 5), ncol = 20)))
long_read <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(rpois(50, 5), ncol = 10)))
SingleCellExperiment::altExp(with_lr, "transcript") <- long_read
SummarizedExperiment::colData(with_lr)$Barcode <- paste0(1:10, "-1")
SummarizedExperiment::colData(without_lr)$Barcode <- paste0(8:27, "-1")
rownames(with_lr) <- as.character(101:110)
rownames(without_lr) <- as.character(103:112)
rownames(long_read) <- as.character(1001:1005)
combined_sce <- FLAMES::combine_sce(sce_with_lr = with_lr, sce_without_lr = without_lr)
#> Missing colnames, using Barcode from colData as colnames
combined_sce
#> class: SingleCellExperiment 
#> dim: 8 27 
#> metadata(0):
#> assays(1): counts
#> rownames(8): 103 104 ... 109 110
#> rowData names(0):
#> colnames: NULL
#> colData names(1): Barcode
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(1): transcript
```
