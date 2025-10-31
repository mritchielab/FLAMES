# Convolution filter for smoothing transcript coverages

Filter out transcripts with sharp drops / rises in coverage, to be used
in `filter_coverage` to remove transcripts with potential misalignments
/ internal priming etc. Filtering is done by convolving the coverage
with a kernal of 1s and -1s (e.g. `c(1, 1, -1, -1)`, where the width of
the 1s and -1s are determined by the `width` parameter), and check if
the maximum absolute value of the convolution is below a threshold. If
the convolution is below the threshold, `TRUE` is returned, otherwise
`FALSE`.

## Usage

``` r
convolution_filter(x, threshold = 0.15, width = 2, trim = 0.05)
```

## Arguments

- x:

  numeric vector of coverage values

- threshold:

  numeric, the threshold for the maximum absolute value of the
  convolution

- width:

  numeric, the width of the 1s and -1s in the kernal. E.g. `width = 2`
  will result in a kernal of `c(1, 1, -1, -1)`

- trim:

  numeric, the proportion of the coverage values to ignore at both ends
  before convolution.

## Value

logical, `TRUE` if the transcript passes the filter, `FALSE` otherwise

## Examples

``` r
# A >30% drop in coverage will fail the filter with threshold = 0.3
convolution_filter(c(1, 1, 1, 0.69, 0.69, 0.69), threshold = 0.3)
#> [1] FALSE
convolution_filter(c(1, 1, 1, 0.71, 0.7, 0.7), threshold = 0.3)
#> [1] TRUE
```
