# Impute missing transcript counts

Impute missing transcript counts using a shared nearest neighbor graph

## Usage

``` r
sc_impute_transcript(combined_sce, dimred = "PCA", ...)
```

## Arguments

- combined_sce:

  A `SingleCellExperiment` object with gene counts and a "transcript"
  altExp slot.

- dimred:

  The name of the reduced dimension to use for building the shared
  nearest neighbor graph.

- ...:

  Additional arguments to pass to
  [`scran::buildSNNGraph`](https://rdrr.io/pkg/scran/man/buildSNNGraph.html).
  E.g. `k = 30`.

## Value

A `SingleCellExperiment` object with imputed logcounts assay in the
"transcript" altExp slot.

## Details

For cells with `NA` values in the "transcript" altExp slot, this
function imputes the missing values from cells with non-missing values.
A shared nearest neighbor graph is built using reduced dimensions from
the `SingleCellExperiment` object, and the imputation is done where the
imputed value for a cell is the weighted sum of the transcript counts of
its neighbors. Imputed values are stored in the "logcounts" assay of the
"transcript" altExp slot. The "counts" assay is used to obtain logcounts
but left unchanged.

## Examples

``` r
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(rpois(50, 5), ncol = 10)))
long_read <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(rpois(40, 5), ncol = 10)))
SingleCellExperiment::altExp(sce, "transcript") <- long_read
SingleCellExperiment::counts(SingleCellExperiment::altExp(sce))[,1:2] <- NA
SingleCellExperiment::counts(SingleCellExperiment::altExp(sce))
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]   NA   NA    1    1    4    2    8    4    4     3
#> [2,]   NA   NA    6    5    7    3    5    6    4     6
#> [3,]   NA   NA    6    3    5    6    9    4    2     8
#> [4,]   NA   NA    4    7    6    7    2    8    5     4
imputed_sce <- sc_impute_transcript(sce, k = 4)
#> Warning: more singular values/vectors requested than available
#> Warning: 'buildSNNGraph' is deprecated.
#> Use 'bluster::makeSNNGraph' instead.
#> See help("Deprecated")
#> Imputing transcript counts ...
SingleCellExperiment::logcounts(SingleCellExperiment::altExp(imputed_sce))
#> 4 x 10 Matrix of class "dgeMatrix"
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#> [1,] 1.828722 1.963206 1.097413 1.144658 2.177193 1.656623 2.898853 2.177193
#> [2,] 2.664177 2.629179 2.970529 2.818582 2.840921 2.080373 2.332410 2.651704
#> [3,] 2.601905 2.545380 2.970529 2.211888 2.433870 2.898853 3.047124 2.177193
#> [4,] 2.709185 2.664554 2.474780 3.244364 2.651704 3.093344 1.386581 3.008174
#>          [,9]    [,10]
#> [1,] 2.624491 1.913744
#> [2,] 2.624491 2.708345
#> [3,] 1.841302 3.067114
#> [4,] 2.898853 2.229734
```
