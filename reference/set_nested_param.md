# Set Nested Configuration Parameter

Helper function to set a nested parameter in a configuration list using
dot notation (e.g., "barcode_parameters.pattern.primer")

## Usage

``` r
set_nested_param(config, param_path, value)
```

## Arguments

- config:

  Configuration list

- param_path:

  Parameter path using dot notation

- value:

  Value to set

## Value

Modified configuration list
