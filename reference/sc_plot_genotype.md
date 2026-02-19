# Plot genotype of single-cell data

Plot the genotype of single-cell data on a reduced dimension plot (e.g.
UMAP).

## Usage

``` r
sc_plot_genotype(
  sce,
  genotype_tb,
  reduced_dim = "UMAP",
  na_cell_col = "grey",
  na_cell_size = 0.1,
  na_cell_alpha = 0.1,
  ...
)
```

## Arguments

- sce:

  SingleCellExperiment: the single-cell experiment object with reduced
  dimensions.

- genotype_tb:

  tibble: the genotype table, output from `sc_genotype`.

- reduced_dim:

  character(1): the name of the reduced dimension to use for plotting.

- na_cell_col:

  character(1): the color of the cells with no genotype.

- na_cell_size:

  numeric(1): the size of the cells with no genotype.

- na_cell_alpha:

  numeric(1): the alpha of the cells with no genotype.

- ...:

  additional arguments passed to `geom_point` for cells with genotype.

## Value

A ggplot2 object with the genotype plotted on the reduced dimension.

## Examples

``` r
ppl <- example_pipeline("SingleCellPipeline") |>
  run_FLAMES()
#> Writing configuration parameters to:  /tmp/RtmpjRxi1E/file95d2784ed24b/config_file_38354.json 
#> Warning: You have set to use oarfish quantification without gene quantification. Oarfish currently does not collapse UMIs, and gene quantification performs UMI collapsing. You may want to set do_gene_quantification to TRUE for more accurate results.
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: FALSE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> ── Running step: barcode_demultiplex @ Thu Feb 19 01:57:04 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/RtmpjRxi1E/file95d2784ed24b/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 8
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
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
#> ── Running step: genome_alignment @ Thu Feb 19 01:57:04 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpjRxi1E/file95d2784ed24b/matched_reads.fastq.gz -> /tmp/RtmpjRxi1E/file95d2784ed24b/align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Thu Feb 19 01:57:04 2026 ─────────────
#> ── Running step: read_realignment @ Thu Feb 19 01:57:05 2026 ───────────────────
#> Checking for fastq file(s) /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpjRxi1E/file95d2784ed24b/matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpjRxi1E/file95d2784ed24b/matched_reads_dedup.fastq.gz
#>  files not found
#> Warning: Oarfish does not support UMI deduplication, you should deduplicate reads before running Oarfish
#> Realigning sample /tmp/RtmpjRxi1E/file95d2784ed24b/matched_reads.fastq.gz -> /tmp/RtmpjRxi1E/file95d2784ed24b/realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> ── Running step: transcript_quantification @ Thu Feb 19 01:57:05 2026 ──────────
sce <- experiment(ppl) |>
 scuttle::logNormCounts() |>
 scater::runPCA() |>
 scater::runUMAP()
#> Warning: more singular values/vectors requested than available
#> using unknown matrix fallback for ' dgTMatrix '
#> Warning: You're computing too large a percentage of total singular values, use a standard svd instead.
snps_tb <- sc_mutations(
  bam_path = ppl@genome_bam,
  seqnames = "chr14",
  positions = 2714
)
#> 01:57:10 Got 1 bam file, parallelizing over each position ...
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#> 
#> 01:57:11 Merging results ...
genotype_tb <- sc_genotype(
  snps_tb, ref = "C", alt = "T", seqname = "chr14", pos = 2714,
  alt_min_count = 2, alt_min_pct = 0.5, ref_min_count = 1, ref_min_pct = 1
)
sc_plot_genotype(
  sce, genotype_tb, na_cell_col = "black",
  na_cell_size = 0.5, na_cell_alpha = 0.7,
  size = 2
)
```
