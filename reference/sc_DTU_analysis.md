# FLAMES Differential Transcript Usage Analysis

Differential transcription usage testing for single cell data, using
`colLabels` as cluster labels.

## Usage

``` r
sc_DTU_analysis(
  sce,
  gene_col = "gene_id",
  min_count = 15,
  threads = 1,
  method = "trascript usage permutation",
  permuations = 1000
)
```

## Arguments

- sce:

  The `SingleCellExperiment` object, with transcript counts in the
  `counts` slot and cluster labels in the `colLabels` slot.

- gene_col:

  The column name in the `rowData` slot of `sce` that contains the gene
  ID / name. Default is `"gene_id"`.

- min_count:

  The minimum total counts for a transcript to be tested.

- threads:

  Number of threads to use for parallel processing.

- method:

  The method to use for testing, listed in `details`.

- permuations:

  Number of permutations for permutation methods.

## Value

a `tibble` containing the following columns:

- p.value:

  \- the raw p-value

- adj.p.value:

  \- multiple testing adjusted p-value

- cluster:

  \- the cluster where DTU was observed

- transcript:

  \- rowname of `sce`, the DTU isoform

- transcript_usage:

  \- the transcript usage of the isoform in the cluster

Additional columns from `method = "trascript usage permutation"`:

- transcript_usage_elsewhere:

  \- transcript usage in other clusters

- usage_difference:

  \- the difference between the two transcript usage

- permuted_var:

  \- the variance of usage difference in the permuted data

Additional columns from `method = "chisq"`:

- X_value:

  \- the test statistic

- df:

  \- the degrees of freedom

- expected_usage:

  \- the expected usage (mean across all clusters)

- usage_difference:

  \- the difference between the observed and expected usage

The table is sorted by P-values.

## Details

Genes with more than 2 isoforms expressing more than `min_count` counts
are selected for testing with one of the following methods:

- trascript usage permutation:

  Transcript usage are taken as the test statistic, cluster labels are
  permuted to generate a null distribution.

- chisq:

  Chi-square test of the transcript count matrix for each gene.

Adjusted P-values were calculated by Benjamini–Hochberg correction.

## Examples

``` r
outdir <- tempfile()
dir.create(outdir)
bc_allow <- file.path(outdir, "bc_allow.tsv")
genome_fa <- file.path(outdir, "rps24.fa")
R.utils::gunzip(
  filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
  destname = bc_allow, remove = FALSE
)
R.utils::gunzip(
  filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
  destname = genome_fa, remove = FALSE
)

sce <- FLAMES::sc_long_pipeline(
  genome_fa = genome_fa,
  fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
  annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  outdir = outdir,
  barcodes_file = bc_allow,
  config_file = create_config(
    outdir,
    pipeline_parameters.demultiplexer = "flexiplex"
  )
)
#> Writing configuration parameters to:  /tmp/RtmpjRxi1E/file95d2734aced8/config_file_38354.json 
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> ── Running step: barcode_demultiplex @ Thu Feb 19 01:56:45 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/RtmpjRxi1E/file95d2734aced8/bc_allow.tsv
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
#> ── Running step: genome_alignment @ Thu Feb 19 01:56:46 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpjRxi1E/file95d2734aced8/matched_reads.fastq.gz -> /tmp/RtmpjRxi1E/file95d2734aced8/align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: gene_quantification @ Thu Feb 19 01:56:46 2026 ────────────────
#> 01:56:46 AM Thu Feb 19 2026 quantify genes 
#> Using BAM(s): /tmp/RtmpjRxi1E/file95d2734aced8/align2genome.bam
#> ── Running step: isoform_identification @ Thu Feb 19 01:56:46 2026 ─────────────
#> ── Running step: read_realignment @ Thu Feb 19 01:56:46 2026 ───────────────────
#> Checking for fastq file(s) /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpjRxi1E/file95d2734aced8/matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpjRxi1E/file95d2734aced8/matched_reads_dedup.fastq.gz
#>  files found
#> Realigning sample /tmp/RtmpjRxi1E/file95d2734aced8/matched_reads_dedup.fastq.gz -> /tmp/RtmpjRxi1E/file95d2734aced8/realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> ── Running step: transcript_quantification @ Thu Feb 19 01:56:47 2026 ──────────
#> Pipeline saved to /tmp/RtmpjRxi1E/file95d2734aced8/pipeline.rds
group_anno <- data.frame(barcode_seq = colnames(sce), groups = SingleCellExperiment::counts(sce)["ENSMUST00000169826.2", ] > 1)
SingleCellExperiment::colLabels(sce) <- group_anno$groups
# DTU with permutation testing:
sc_DTU_analysis(sce, min_count = 1, method = "trascript usage permutation")
#> Filtering for genes with at least 2 isforms expressing more than 1 counts ...
#>  1 gene(s), 8 transcript(s) left.
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#> 
#> # A tibble: 16 × 9
#>    transcript   cluster transcript_usage transcript_usage_els…¹ usage_difference
#>    <chr>        <fct>              <dbl>                  <dbl>            <dbl>
#>  1 ENSMUSG0000… FALSE             0.0909                 0.0261         0.0648  
#>  2 ENSMUSG0000… TRUE              0.0261                 0.0909        -0.0648  
#>  3 ENSMUST0000… FALSE             0.608                  0.755         -0.147   
#>  4 ENSMUST0000… TRUE              0.755                  0.608          0.147   
#>  5 ENSMUSG0000… FALSE             0.0568                 0.0201         0.0367  
#>  6 ENSMUSG0000… TRUE              0.0201                 0.0568        -0.0367  
#>  7 ENSMUSG0000… FALSE             0.119                  0.0909         0.0284  
#>  8 ENSMUSG0000… FALSE             0.0341                 0.0234         0.0106  
#>  9 ENSMUSG0000… FALSE             0.0568                 0.0446         0.0122  
#> 10 ENSMUSG0000… TRUE              0.0909                 0.119         -0.0284  
#> 11 ENSMUSG0000… TRUE              0.0234                 0.0341        -0.0106  
#> 12 ENSMUSG0000… TRUE              0.0446                 0.0568        -0.0122  
#> 13 ENSMUSG0000… FALSE             0.0227                 0.0288        -0.00607 
#> 14 ENSMUSG0000… TRUE              0.0288                 0.0227         0.00607 
#> 15 ENSMUST0000… FALSE             0.0114                 0.0111         0.000224
#> 16 ENSMUST0000… TRUE              0.0111                 0.0114        -0.000224
#> # ℹ abbreviated name: ¹​transcript_usage_elsewhere
#> # ℹ 4 more variables: p.value <dbl>, permuted_var <dbl>, adj.p.value <dbl>,
#> #   gene_id <chr>
# now try with chisq:
sc_DTU_analysis(sce, min_count = 1, method = "chisq")
#> Filtering for genes with at least 2 detected isforms ...
#>  8 isoform(s) left.
#> Aggregating counts by cluster labels ...
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#> 
#> # A tibble: 1 × 10
#>   gene  X_value    df transcript cluster p.value expected_usage transcript_usage
#>   <chr>   <dbl> <int> <chr>      <chr>     <dbl>          <dbl>            <dbl>
#> 1 ENSM…    12.8     7 ENSMUSG00… FALSE    0.0769         0.0421           0.0909
#> # ℹ 2 more variables: usage_difference <dbl>, adj.p.value <dbl>
```
