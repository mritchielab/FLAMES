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
#> ℹ Writing configuration to: /tmp/RtmpnC89xy/filebc1526028c7d/config_file_48149.json
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> ── Running step: barcode_demultiplex @ Thu Aug  6 04:49:59 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/RtmpnC89xy/filebc1526028c7d/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 8
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
#> ── Running step: genome_alignment @ Thu Aug  6 04:50:00 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpnC89xy/filebc1526028c7d/matched_reads.fastq.gz -> /tmp/RtmpnC89xy/filebc1526028c7d/align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: gene_quantification @ Thu Aug  6 04:50:00 2026 ────────────────
#> 04:50:00 AM Thu Aug 06 2026 quantify genes 
#> Using BAM(s): /tmp/RtmpnC89xy/filebc1526028c7d/align2genome.bam
#> ── Running step: isoform_identification @ Thu Aug  6 04:50:00 2026 ─────────────
#> ── Running step: read_realignment @ Thu Aug  6 04:50:00 2026 ───────────────────
#> Checking for fastq file(s) /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpnC89xy/filebc1526028c7d/matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpnC89xy/filebc1526028c7d/matched_reads_dedup.fastq.gz
#>  files found
#> Realigning sample /tmp/RtmpnC89xy/filebc1526028c7d/matched_reads_dedup.fastq.gz -> /tmp/RtmpnC89xy/filebc1526028c7d/realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> ── Running step: transcript_quantification @ Thu Aug  6 04:50:01 2026 ──────────
sce <- experiment(ppl) |>
 scrapper::normalizeRnaCounts.se() |>
 scater::runPCA() |>
 scater::runUMAP()
#> using unknown matrix fallback for 'dgTMatrix'
#> Warning: more singular values/vectors requested than available
#> using unknown matrix fallback for 'dgTMatrix'
snps_tb <- sc_mutations(
  bam_path = ppl@genome_bam,
  seqnames = "chr14",
  positions = 2714
)
#> 04:50:06 Got 1 bam file, parallelizing over each position ...
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#> 
#> 04:50:07 Merging results ...
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
