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
#> [1,]   NA   NA    2    8    3    4    5    5    5     8
#> [2,]   NA   NA    7    5    3    7    4    3    8     5
#> [3,]   NA   NA    3    3    5    7    3    9    6     8
#> [4,]   NA   NA    2    3    6    7    4    1    4     5
imputed_sce <- sc_impute_transcript(sce, k = 4)
#> Warning: more singular values/vectors requested than available
#> Warning: You're computing too large a percentage of total singular values, use a standard svd instead.
#> Imputing transcript counts ...
SingleCellExperiment::logcounts(SingleCellExperiment::altExp(imputed_sce))
#> 4 x 10 Matrix of class "dgeMatrix"
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#> [1,] 2.513033 2.482428 1.934112 3.219678 2.165203 2.056584 2.842350 2.697354
#> [2,] 2.616786 2.568155 3.442943 2.631656 2.165203 2.707083 2.569856 2.101538
#> [3,] 2.589690 2.699316 2.387402 2.042091 2.767406 2.707083 2.233620 3.442943
#> [4,] 2.269120 2.162505 1.934112 2.042091 2.994686 2.707083 2.569856 1.068480
#>          [,9]    [,10]
#> [1,] 2.404216 2.823122
#> [2,] 2.976284 2.262456
#> [3,] 2.621096 2.823122
#> [4,] 2.148863 2.262456
```
