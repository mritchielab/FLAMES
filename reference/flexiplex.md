# Rcpp port of flexiplex

demultiplex reads with flexiplex, for detailed description, see
documentation for the original flexiplex:
https://davidsongroup.github.io/flexiplex

## Usage

``` r
flexiplex(
  reads_in,
  barcodes_file,
  bc_as_readid,
  max_bc_editdistance,
  max_flank_editdistance,
  pattern,
  reads_out,
  stats_out,
  bc_out,
  reverseCompliment,
  n_threads
)
```

## Arguments

- reads_in:

  Input FASTQ or FASTA file

- barcodes_file:

  barcode allow-list file

- bc_as_readid:

  bool, whether to add the demultiplexed barcode to the read ID field

- max_bc_editdistance:

  max edit distance for barcode '

- max_flank_editdistance:

  max edit distance for the flanking sequences '

- pattern:

  StringVector defining the barcode structure, see \[find_barcode\]

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
