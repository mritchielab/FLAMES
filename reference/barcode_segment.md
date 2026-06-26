# Create a Flexiplex barcode segment

Creates a `FlexiplexSegment` object describing one component of a read's
barcode structure. A list of segments is passed to
[`find_barcode`](https://mritchielab.github.io/FLAMES/reference/find_barcode.md)
(via the `segments` argument) to define how reads are parsed from 5' to
3'.

## Usage

``` r
barcode_segment(
  type = "FIXED",
  pattern,
  name,
  bc_list = NA_character_,
  group = NA_character_,
  buffer_size = 2,
  max_edit_distance = 2
)
```

## Arguments

- type:

  Segment type: one of `"FIXED"`, `"MATCHED"`, `"RANDOM"`, or
  `"MATCHED_SPLIT"`.

- pattern:

  The nucleotide pattern for this segment. For `FIXED`: the known
  sequence (e.g. `"CTACACGACGCTCTTCCGATCT"`). For
  `MATCHED`/`MATCHED_SPLIT`: an N-repeat matching the expected barcode
  length (e.g. `"NNNNNNNNNNNNNNNN"` for a 16-nt barcode). For `RANDOM`:
  an N-repeat matching the UMI length (e.g. `"NNNNNNNNNNNN"` for a 12-nt
  UMI). Other IUPAC ambigious codes such as "Y" "R" are also supported,
  e.g. `"NNNYNNNRNNN"`, `"CTACACGACGCTCTTCCGATCTNNN"`

- name:

  Label for this segment in output files (e.g. `"CB"` for cell barcode,
  `"UB"` for UMI, `"primer"`). Defaults to `"FIXED_SEGMENT"` for `FIXED`
  segments.

- bc_list:

  For `MATCHED` segments: a key used to look up the barcode allow-list
  file from the `barcodes_files` argument of
  [`find_barcode`](https://mritchielab.github.io/FLAMES/reference/find_barcode.md).
  The key must match a name in `barcodes_files`, or `barcodes_files` can
  be a single unnamed path (in which case `bc_list` may be omitted).

- group:

  For `MATCHED_SPLIT` segments: the name of the
  [`barcode_group`](https://mritchielab.github.io/FLAMES/reference/barcode_group.md)
  this segment belongs to. Must match the `name` of a
  [`barcode_group`](https://mritchielab.github.io/FLAMES/reference/barcode_group.md)
  object passed to
  [`find_barcode`](https://mritchielab.github.io/FLAMES/reference/find_barcode.md).

- buffer_size:

  Non-negative integer. Extra nucleotides searched on each side of the
  expected segment position to accommodate small insertions/deletions.
  Only applies to `MATCHED` and `MATCHED_SPLIT` segments. Default: `2`.

- max_edit_distance:

  Non-negative integer. Maximum edit distance allowed when matching a
  read sequence against the barcode allow-list. Only applies to
  `MATCHED` and `MATCHED_SPLIT` segments. Default: `2`.

## Value

A `FlexiplexSegment` object for use in
[`find_barcode`](https://mritchielab.github.io/FLAMES/reference/find_barcode.md).

## Details

Four segment types are supported:

- `"FIXED"`:

  A known, constant flanking sequence (e.g. a sequencing primer or
  poly-T tail). Used as an alignment anchor; i.e. no barcode allow-list
  is associated

- `"MATCHED"`:

  A variable sequence matched against a barcode allow-list (e.g. a cell
  barcode). The allow-list file is resolved via the `barcodes_files`
  argument of
  [`find_barcode`](https://mritchielab.github.io/FLAMES/reference/find_barcode.md).
  `bc_list` must be supplied (or `barcodes_files` must be a single
  unnamed path).

- `"RANDOM"`:

  A random sequence of fixed length captured verbatim without matching
  (e.g. a UMI). No allow-list is needed.

- `"MATCHED_SPLIT"`:

  Like `"MATCHED"`, but participates in a
  [`barcode_group`](https://mritchielab.github.io/FLAMES/reference/barcode_group.md)
  for multi-segment barcode matching, where the all MATCHED_SPLIT
  segments of the same group are concatenated and then matched to the
  allow-list barcodes. `group` must name a corresponding
  [`barcode_group`](https://mritchielab.github.io/FLAMES/reference/barcode_group.md).

## See also

[`barcode_group`](https://mritchielab.github.io/FLAMES/reference/barcode_group.md),
[`find_barcode`](https://mritchielab.github.io/FLAMES/reference/find_barcode.md)

## Examples

``` r
# A typical 10x Genomics 3' v3 barcode structure:
segments <- list(
  barcode_segment(type = "FIXED",   pattern = "CTACACGACGCTCTTCCGATCT", name = "primer"),
  barcode_segment(type = "MATCHED", pattern = "NNNNNNNNNNNNNNNN",        name = "CB",
                  bc_list = "CB", buffer_size = 5, max_edit_distance = 2),
  barcode_segment(type = "RANDOM",  pattern = "NNNNNNNNNNNN",            name = "UB"),
  barcode_segment(type = "FIXED",   pattern = "TTTTTTTTT",               name = "polyT")
)
```
