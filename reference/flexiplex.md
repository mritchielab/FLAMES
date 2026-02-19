# Rcpp port of flexiplex

demultiplex reads with flexiplex, for detailed description, see
documentation for the original flexiplex:
https://davidsongroup.github.io/flexiplex

## Usage

``` r
flexiplex(
  r_segments,
  r_barcode_groups,
  max_flank_editdistance,
  reads_in,
  reads_out,
  stats_out,
  bc_out,
  reverseCompliment,
  n_threads
)
```

## Arguments

- r_segments:

  List defining the barcode structure

- r_barcode_groups:

  List defining barcode groups

- max_flank_editdistance:

  int, maximum edit distance for matching flanking sequences

- reads_in:

  Input FASTQ or FASTA file

- reads_out:

  output file for demultiplexed reads

- stats_out:

  output file for demultiplexed stats

- bc_out:

  WIP

- reverseCompliment:

  bool, whether to reverse complement the reads after demultiplexing

- n_threads:

  number of threads to be used during demultiplexing

## Value

integer return value. 0 represents normal return.
