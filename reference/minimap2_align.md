# Minimap2 Align to Genome

Uses minimap2 to align sequences agains a reference databse. Uses
options '-ax splice -t 12 -k14 –secondary=no `fa_file` `fq_in`'

## Usage

``` r
minimap2_align(
  fq_in,
  fa_file,
  config,
  outfile,
  minimap2_args,
  sort_by,
  minimap2,
  samtools,
  threads = 1,
  tmpdir
)
```

## Arguments

- fq_in:

  File path to the fastq file used as a query sequence file

- fa_file:

  Path to the fasta file used as a reference database for alignment

- config:

  Parsed list of FLAMES config file

- outfile:

  Path to the output file

- minimap2_args:

  Arguments to pass to minimap2, see minimap2 documentation for details.

- sort_by:

  Column to sort the bam file by, see `samtools sort` for details

- minimap2:

  Path to minimap2 binary

- samtools:

  path to the samtools binary.

- threads:

  Integer, threads for minimap2 to use, see minimap2 documentation for
  details,

- tmpdir:

  Temporary directory to use for intermediate files. FLAMES will try to
  detect cores if this parameter is not provided.

## Value

a `data.frame` summarising the reads aligned
