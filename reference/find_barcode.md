# Match Cell Barcodes

demultiplex reads with flexiplex

## Usage

``` r
find_barcode(
  fastq,
  barcodes_file,
  max_bc_editdistance = 2,
  max_flank_editdistance = 8,
  reads_out,
  stats_out,
  threads = 1,
  pattern = c(primer = "CTACACGACGCTCTTCCGATCT", BC = paste0(rep("N", 16), collapse =
    ""), UMI = paste0(rep("N", 12), collapse = ""), polyT = paste0(rep("T", 9), collapse
    = "")),
  TSO_seq = "",
  TSO_prime = 3,
  strand = "+",
  cutadapt_minimum_length = 1,
  full_length_only = FALSE
)
```

## Arguments

- fastq:

  A path to a FASTQ file or a directory containing FASTQ files.

- barcodes_file:

  path to file containing barcode allow-list, with one barcode in each
  line

- max_bc_editdistance:

  max edit distances for the barcode sequence

- max_flank_editdistance:

  max edit distances for the flanking sequences (primer and polyT)

- reads_out:

  path to output FASTQ file

- stats_out:

  path of output stats file

- threads:

  number of threads to be used

- pattern:

  named character vector defining the barcode pattern

- TSO_seq:

  TSO sequence to be trimmed

- TSO_prime:

  either 3 (when `TSO_seq` is on 3' the end) or 5 (on 5' end)

- strand:

  strand of the barcode pattern, either '+' or '-' (read will be reverse
  complemented after barcode matching if '-')

- cutadapt_minimum_length:

  minimum read length after TSO trimming (cutadapt's –minimum-length)

- full_length_only:

  boolean, when TSO sequence is provided, whether reads without TSO are
  to be discarded

## Value

a list containing: `reads_tb` (tibble of read demultiplexed information)
and `input`, `output`, `read1_with_adapter` from cutadapt report (if TSO
trimming is performed)

## Details

This function demultiplexes reads by searching for flanking sequences
(adaptors) around the barcode sequence, and then matching against
allowed barcodes.

## Examples

``` r
outdir <- tempfile()
dir.create(outdir)
bc_allow <- file.path(outdir, "bc_allow.tsv")
R.utils::gunzip(
  filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
  destname = bc_allow, remove = FALSE
)
find_barcode(
  fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
  stats_out = file.path(outdir, "bc_stat.tsv.gz"),
  reads_out = file.path(outdir, "demultiplexed.fastq.gz"),
  barcodes_file = bc_allow, 
  TSO_seq = "AAGCAGTGGTATCAACGCAGAGTACATGGG", TSO_prime = 5,
  strand = '-', cutadapt_minimum_length = 10, full_length_only = TRUE
)
#> FLEXIPLEX 0.96.2
#> Setting max barcode edit distance to 2
#> Setting max flanking sequence edit distance to 8
#> Setting read IDs to be  replaced
#> Setting number of threads to 1
#> Search pattern: 
#> primer: CTACACGACGCTCTTCCGATCT
#> BC: NNNNNNNNNNNNNNNN
#> UMI: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Setting known barcodes from /tmp/RtmpaNKtnp/file80d01712335e/bc_allow.tsv
#> Number of known barcodes: 143
#> Processing file: /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> Searching for barcodes...
#> Number of reads processed: 393
#> Number of reads where at least one barcode was found: 368
#> Number of reads with exactly one barcode match: 364
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
#>                                       total reads 
#>                                               393 
#>                               demultiplexed reads 
#>                                               368 
#>                                single match reads 
#>                                               364 
#> single strand single barcode multi-matching reads 
#>                                                 0 
#>              single strand multiple barcode reads 
#>                                                 3 
#>                 both strands single barcode reads 
#>                                                 0 
#>               both strands multiple barcode reads 
#>                                                 1 
#> 
#> $stats_out
#> [1] "/tmp/RtmpaNKtnp/file80d01712335e/bc_stat.tsv.gz"
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
#>  [4] "/tmp/RtmpaNKtnp/file80d01712335e/demultiplexed.fastq.gz"          
#>  [5] "/tmp/RtmpaNKtnp/file80d01712335e/untrimmed_demultiplexed.fastq.gz"
#>  [6] "--json"                                                           
#>  [7] "/tmp/RtmpaNKtnp/file80d01712335e/file80d015f4419c.json"           
#>  [8] "--untrimmed-output"                                               
#>  [9] "/tmp/RtmpaNKtnp/file80d01712335e/noTSO_demultiplexed.fastq.gz"    
#> [10] "--minimum-length"                                                 
#> [11] "10"                                                               
#> 
#> $cutadapt$cores
#> [1] 1
#> 
#> $cutadapt$input
#> $cutadapt$input$path1
#> [1] "/tmp/RtmpaNKtnp/file80d01712335e/untrimmed_demultiplexed.fastq.gz"
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
#> [1] 248786
#> 
#> $cutadapt$basepair_counts$input_read1
#> [1] 248786
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
#> [1] 154283
#> 
#> $cutadapt$basepair_counts$output_read1
#> [1] 154283
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
