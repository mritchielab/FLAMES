# Parse FLAMES' GFF output

Parse FLAMES' GFF ouputs into a Genomic Ranges List

## Usage

``` r
get_GRangesList(
  file,
  feature.type = c("exon", "utr"),
  drop.cols = c("type", "exon_number", "exon_id", "level", "Parent")
)
```

## Arguments

- file:

  the GFF file to parse

- feature.type:

  The type of features to extract from the GFF file. Default is
  `c("exon", "utr")`.

- drop.cols:

  Columns to drop from the metadata. Default is
  `c("type", "exon_number", "exon_id", "level")`, which are
  exon-specific metadata that may not be relevant when keeping just the
  first row (exon).

## Value

A list containing a `GRangesList` of isoforms and a `DataFrame`, which
have the same number of rows as the number of unique transcript IDs in
the GFF file.
