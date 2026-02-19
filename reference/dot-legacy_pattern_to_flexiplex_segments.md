# Convert legacy pattern to Flexiplex segments

Convert legacy pattern to Flexiplex segments

## Usage

``` r
.legacy_pattern_to_flexiplex_segments(
  pattern,
  barcodes_file,
  max_bc_editdistance,
  buffer_size
)
```

## Arguments

- pattern:

  named character vector defining the barcode pattern

- barcodes_file:

  path to file containing barcode allow-list

- max_bc_editdistance:

  max edit distances for the barcode sequence

- buffer_size:

  buffer size for barcode matching

## Value

a list of FlexiplexSegment objects
