# Load Configurations

Loads a configuration file and fills in missing values with defaults
from the package's default configuration.

## Usage

``` r
load_config(config_file, type = "sc_3end")
```

## Arguments

- config_file:

  Path to the configuration JSON file

- type:

  Config type to use for defaults ("sc_3end" or "SIRV")

## Value

A complete configuration list with all parameters filled
