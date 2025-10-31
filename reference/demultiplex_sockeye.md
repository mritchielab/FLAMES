# Demultiplex reads using Sockeye outputs

Demultiplex reads using the `cell_umi_gene.tsv` file from Sockeye.

## Usage

``` r
demultiplex_sockeye(fastq_dir, sockeye_tsv, out_fq)
```

## Arguments

- fastq_dir:

  The folder containing FASTQ files from Sockeye's output under
  `ingest/chunked_fastqs`.

- sockeye_tsv:

  The `cell_umi_gene.tsv` file from Sockeye.

- out_fq:

  The output FASTQ file.

## Value

returns NULL
