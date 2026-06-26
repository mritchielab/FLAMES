# Match Cell Barcodes

demultiplex reads with flexiplex

## Usage

``` r
find_barcode(
  fastq,
  segments,
  barcode_groups,
  barcodes_files,
  max_flank_editdistance = 8,
  reads_out,
  stats_out,
  threads = 1,
  TSO_seq = "",
  TSO_prime = 3,
  strand = "+",
  cutadapt_minimum_length = 1,
  full_length_only = FALSE,
  pattern = c(primer = "CTACACGACGCTCTTCCGATCT", BC = paste0(rep("N", 16), collapse =
    ""), UMI = paste0(rep("N", 12), collapse = ""), polyT = paste0(rep("T", 9), collapse
    = "")),
  max_bc_editdistance = 2
)
```

## Arguments

- fastq:

  A path to a FASTQ file or a directory containing FASTQ files.

- segments:

  A list of
  [`barcode_segment`](https://mritchielab.github.io/FLAMES/reference/barcode_segment.md)
  (`FlexiplexSegment`) objects describing the read structure. May also
  be a data frame with columns `type`, `pattern`, `name`, and optionally
  `bc_list_name`, `group`, `buffer_size`, and `max_edit_distance` (as
  produced by `jsonlite` when reading a config file). When omitted, the
  legacy `pattern` argument is used instead.

- barcode_groups:

  A list of
  [`barcode_group`](https://mritchielab.github.io/FLAMES/reference/barcode_group.md)
  (`FlexiplexGroup`) objects for multi-segment barcode matching.
  Required only when `MATCHED_SPLIT` segments are present; pass an empty
  list [`list()`](https://rdrr.io/r/base/list.html) otherwise.

- barcodes_files:

  Path(s) to barcode allow-list file(s), one barcode per line. Can be:

  - A single unnamed character string — used for all `MATCHED` segments
    and barcode groups.

  - A named character vector — names must match the `bc_list` keys set
    in
    [`barcode_segment`](https://mritchielab.github.io/FLAMES/reference/barcode_segment.md)
    and the `bc_list_name` keys set in
    [`barcode_group`](https://mritchielab.github.io/FLAMES/reference/barcode_group.md).

- max_flank_editdistance:

  Maximum edit distance for matching the fixed flanking sequences (e.g.
  primer, poly-T tail). Default: `8`.

- reads_out:

  Path to output FASTQ file containing demultiplexed reads with barcode
  information added to the read header.

- stats_out:

  Path to output TSV (optionally gzip-compressed) file with per-read
  demultiplex results.

- threads:

  Number of threads to use. Default: `1`.

- TSO_seq:

  TSO (template-switching oligo) sequence to trim from reads using
  cutadapt. Set to `""` (default) to skip TSO trimming.

- TSO_prime:

  Either `5` or `3`, indicating whether `TSO_seq` is located at the 5'
  or 3' end of the read after barcode demultiplexing.

- strand:

  Strand of the barcode pattern, either `"+"` or `"-"`. When `"-"`,
  reads are reverse-complemented after barcode matching so that the
  transcript sequence is in the sense direction. Default: `"+"`.

- cutadapt_minimum_length:

  Minimum read length (in nucleotides) to retain after TSO trimming
  (passed to cutadapt's `--minimum-length`). Default: `1`.

- full_length_only:

  Logical. When `TSO_seq` is provided, whether to discard reads that do
  not contain the TSO (i.e. keep only full-length reads). Default:
  `FALSE`.

- pattern:

  **Deprecated.** Named character vector defining the barcode structure
  (legacy interface). Use `segments` instead. Entries named `"BC"` are
  treated as `MATCHED` segments, `"UMI"` as `RANDOM`, and all others as
  `FIXED`.

- max_bc_editdistance:

  **Deprecated.** Maximum edit distance for barcode matching when using
  the legacy `pattern` argument.

## Value

A list containing:

- `read_counts`: a named integer vector with `total reads`,
  `demultiplexed reads`, and `single match reads`.

- `stats_out`: the path to the demultiplex stats file.

- `cutadapt`: (only when TSO trimming is performed) the parsed cutadapt
  JSON report.

## Details

This function demultiplexes reads by searching for flanking sequences
(adaptors) around the barcode sequence, and then matching against an
allowed barcode list.

The read structure is described by a list of
[`barcode_segment`](https://mritchielab.github.io/FLAMES/reference/barcode_segment.md)
objects passed to `segments`. Each segment describes one component of
the read (e.g. a fixed primer, a cell barcode, a UMI, a poly-T tail).
Use
[`barcode_segment`](https://mritchielab.github.io/FLAMES/reference/barcode_segment.md)
to construct segments and, for combinatorial multi-segment barcodes,
[`barcode_group`](https://mritchielab.github.io/FLAMES/reference/barcode_group.md)
to define groups.

For backward compatibility, the legacy `pattern` argument (a named
character vector) is still accepted when `segments` is not supplied.

## See also

[`barcode_segment`](https://mritchielab.github.io/FLAMES/reference/barcode_segment.md),
[`barcode_group`](https://mritchielab.github.io/FLAMES/reference/barcode_group.md),
[`plot_demultiplex_raw`](https://mritchielab.github.io/FLAMES/reference/plot_demultiplex_raw.md)

## Examples

``` r
outdir <- tempfile()
dir.create(outdir)
bc_allow <- file.path(outdir, "bc_allow.tsv")
R.utils::gunzip(
  filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
  destname = bc_allow, remove = FALSE
)

# Modern interface: define segments explicitly
find_barcode(
  fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
  segments = list(
    barcode_segment("FIXED",   "CTACACGACGCTCTTCCGATCT", "primer"),
    barcode_segment("MATCHED", "NNNNNNNNNNNNNNNN", "CB",
                    bc_list = "CB", buffer_size = 5, max_edit_distance = 2),
    barcode_segment("RANDOM",  "NNNNNNNNNNNN", "UB"),
    barcode_segment("FIXED",   "TTTTTTTTT", "polyT")
  ),
  barcode_groups = list(),
  barcodes_files = c(CB = bc_allow),
  stats_out = file.path(outdir, "bc_stat.tsv.gz"),
  reads_out = file.path(outdir, "demultiplexed.fastq.gz"),
  TSO_seq = "AAGCAGTGGTATCAACGCAGAGTACATGGG", TSO_prime = 5,
  strand = '-', cutadapt_minimum_length = 10, full_length_only = TRUE
)
#> Loading known barcodes from /tmp/Rtmp4nGYdi/filebc445d9d090b/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> CB:Z: tag field: CB
#> Processing file: /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> Searching for barcodes...
#> Number of reads processed: 393
#> Number of reads where at least one barcode was found: 368
#> Number of chimera reads: 1
#> All done!
#> Reads    Barcodes
#> 10   2
#> 9    2
#> 8    5
#> 7    4
#> 6    3
#> 5    7
#> 4    14
#> 3    14
#> 2    29
#> 1    57
#> $read_counts
#>         total reads demultiplexed reads  single match reads       chimera reads 
#>                 393                 368                 364                   1 
#> 
#> $stats_out
#> [1] "/tmp/Rtmp4nGYdi/filebc445d9d090b/bc_stat.tsv.gz"
#> 
#> $cutadapt
#> $cutadapt$tag
#> [1] "Cutadapt report"
#> 
#> $cutadapt$schema_version
#> [1] 0 3
#> 
#> $cutadapt$cutadapt_version
#> [1] "4.9"
#> 
#> $cutadapt$python_version
#> [1] "3.11.9"
#> 
#> $cutadapt$command_line_arguments
#>  [1] "-g"                                                               
#>  [2] "AAGCAGTGGTATCAACGCAGAGTACATGGG"                                   
#>  [3] "-o"                                                               
#>  [4] "/tmp/Rtmp4nGYdi/filebc445d9d090b/demultiplexed.fastq.gz"          
#>  [5] "/tmp/Rtmp4nGYdi/filebc445d9d090b/untrimmed_demultiplexed.fastq.gz"
#>  [6] "--json"                                                           
#>  [7] "/tmp/Rtmp4nGYdi/filebc445d9d090b/filebc444978f443.json"           
#>  [8] "--untrimmed-output"                                               
#>  [9] "/tmp/Rtmp4nGYdi/filebc445d9d090b/noTSO_demultiplexed.fastq.gz"    
#> [10] "--minimum-length"                                                 
#> [11] "10"                                                               
#> 
#> $cutadapt$cores
#> [1] 1
#> 
#> $cutadapt$input
#> $cutadapt$input$path1
#> [1] "/tmp/Rtmp4nGYdi/filebc445d9d090b/untrimmed_demultiplexed.fastq.gz"
#> 
#> $cutadapt$input$path2
#> NULL
#> 
#> $cutadapt$input$paired
#> [1] FALSE
#> 
#> 
#> $cutadapt$read_counts
#> $cutadapt$read_counts$input
#> [1] 372
#> 
#> $cutadapt$read_counts$filtered
#> $cutadapt$read_counts$filtered$too_short
#> [1] 0
#> 
#> $cutadapt$read_counts$filtered$too_long
#> NULL
#> 
#> $cutadapt$read_counts$filtered$too_many_n
#> NULL
#> 
#> $cutadapt$read_counts$filtered$too_many_expected_errors
#> NULL
#> 
#> $cutadapt$read_counts$filtered$casava_filtered
#> NULL
#> 
#> $cutadapt$read_counts$filtered$discard_trimmed
#> NULL
#> 
#> $cutadapt$read_counts$filtered$discard_untrimmed
#> NULL
#> 
#> 
#> $cutadapt$read_counts$output
#> [1] 248
#> 
#> $cutadapt$read_counts$reverse_complemented
#> NULL
#> 
#> $cutadapt$read_counts$read1_with_adapter
#> [1] 248
#> 
#> $cutadapt$read_counts$read2_with_adapter
#> NULL
#> 
#> 
#> $cutadapt$basepair_counts
#> $cutadapt$basepair_counts$input
#> [1] 251665
#> 
#> $cutadapt$basepair_counts$input_read1
#> [1] 251665
#> 
#> $cutadapt$basepair_counts$input_read2
#> NULL
#> 
#> $cutadapt$basepair_counts$quality_trimmed
#> NULL
#> 
#> $cutadapt$basepair_counts$quality_trimmed_read1
#> NULL
#> 
#> $cutadapt$basepair_counts$quality_trimmed_read2
#> NULL
#> 
#> $cutadapt$basepair_counts$poly_a_trimmed
#> NULL
#> 
#> $cutadapt$basepair_counts$poly_a_trimmed_read1
#> NULL
#> 
#> $cutadapt$basepair_counts$poly_a_trimmed_read2
#> NULL
#> 
#> $cutadapt$basepair_counts$output
#> [1] 156234
#> 
#> $cutadapt$basepair_counts$output_read1
#> [1] 156234
#> 
#> $cutadapt$basepair_counts$output_read2
#> NULL
#> 
#> 
#> $cutadapt$adapters_read1
#>   name total_matches on_reverse_complement linked five_prime_end.type
#> 1    1           248                    NA  FALSE  regular_five_prime
#>          five_prime_end.sequence five_prime_end.error_rate
#> 1 AAGCAGTGGTATCAACGCAGAGTACATGGG                       0.1
#>   five_prime_end.indels five_prime_end.error_lengths five_prime_end.matches
#> 1                  TRUE                 9, 19, 2....                    248
#>   five_prime_end.adjacent_bases five_prime_end.dominant_adjacent_base
#> 1                            NA                                    NA
#>   five_prime_end.trimmed_lengths three_prime_end
#> 1                   c(3, 16,....              NA
#> 
#> $cutadapt$adapters_read2
#> NULL
#> 
#> $cutadapt$poly_a_trimmed_read1
#> NULL
#> 
#> $cutadapt$poly_a_trimmed_read2
#> NULL
#> 
#> 
```
