# Recursively Merge Configuration Lists

Internal function to recursively merge configuration lists, filling
missing values from defaults while preserving user values

## Usage

``` r
merge_configs_recursive(default_config, user_config)
```

## Arguments

- default_config:

  Default configuration list

- user_config:

  User configuration list

## Value

Merged configuration list

## Note

Special case: when user_config contains barcode_parameters.pattern as a
list, the entire pattern list is preserved as-is without merging with
defaults to maintain user-specified order and structure.
