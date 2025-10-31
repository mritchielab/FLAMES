# Find path to a binary Wrapper for Sys.which to find path to a binary

This function is a wrapper for
[`base::Sys.which`](https://rdrr.io/r/base/Sys.which.html) to find the
path to a command. It also searches within the `FLAMES` basilisk conda
environment. This function also replaces "" with `NA` in the output of
[`base::Sys.which`](https://rdrr.io/r/base/Sys.which.html) to make it
easier to check if the binary is found.

## Usage

``` r
find_bin(command)
```

## Arguments

- command:

  character, the command to search for

## Value

character, the path to the command or `NA`

## Examples

``` r
find_bin("minimap2")
#>                                 minimap2 
#> "/__w/_temp/Library/FLAMES/bin/minimap2" 
```
