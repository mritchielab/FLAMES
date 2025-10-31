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
#> [1,]   NA   NA    5    1    3    8    4    3    7     7
#> [2,]   NA   NA    6    3    4    1    6    4    5     2
#> [3,]   NA   NA    3   12    5    4    3    6    9     4
#> [4,]   NA   NA    7    8    8    5    3    6    9     6
imputed_sce <- sc_impute_transcript(sce, k = 4)
#> Warning: more singular values/vectors requested than available
#> Warning: You're computing too large a percentage of total singular values, use a standard svd instead.
#> Imputing transcript counts ...
SingleCellExperiment::logcounts(SingleCellExperiment::altExp(imputed_sce))
#> 4 x 10 Matrix of class "dgeMatrix"
#>          [,1]     [,2]     [,3]      [,4]     [,5]     [,6]     [,7]     [,8]
#> [1,] 2.394844 2.474454 2.577788 0.9028775 2.046578 3.361456 2.636625 2.103012
#> [2,] 2.145317 2.157884 2.799975 1.8517490 2.371559 1.110846 3.142107 2.431553
#> [3,] 2.700397 2.596905 1.993545 3.5156998 2.636625 2.495411 2.296916 2.924500
#> [4,] 2.861925 2.910323 2.992466 2.9924663 3.224966 2.765240 2.296916 2.924500
#>          [,9]    [,10]
#> [1,] 2.553565 3.119487
#> [2,] 2.163230 1.676885
#> [3,] 2.860466 2.431553
#> [4,] 2.860466 2.924500
```
