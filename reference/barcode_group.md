# Create a Flexiplex barcode group

Creates a `FlexiplexGroup` object that links together a set of
`MATCHED_SPLIT`
[`barcode_segment`](https://mritchielab.github.io/FLAMES/reference/barcode_segment.md)s
for concatenatiing segments before barcode matching. A barcode group
lets flexiplex jointly match multiple segment sequences against a shared
allow-list (e.g. the concatenation of two barcode halves).

## Usage

``` r
barcode_group(name, bc_list_name, max_edit_distance = 2)
```

## Arguments

- name:

  The group identifier. Must match the `group` argument of every
  [`barcode_segment`](https://mritchielab.github.io/FLAMES/reference/barcode_segment.md)
  that belongs to this group.

- bc_list_name:

  A key used to look up the allow-list file for this group from the
  `barcodes_files` argument of
  [`find_barcode`](https://mritchielab.github.io/FLAMES/reference/find_barcode.md).
  The key must match a name in `barcodes_files`.

- max_edit_distance:

  Non-negative integer. Maximum total edit distance allowed when
  matching the concatenated segment sequences against the group
  allow-list. Default: `2`.

## Value

A `FlexiplexGroup` object for use in
[`find_barcode`](https://mritchielab.github.io/FLAMES/reference/find_barcode.md).

## See also

[`barcode_segment`](https://mritchielab.github.io/FLAMES/reference/barcode_segment.md),
[`find_barcode`](https://mritchielab.github.io/FLAMES/reference/find_barcode.md)
