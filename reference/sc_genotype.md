# Genotype a single-cell mutation

A simplistic function to genotype a single-cell mutation at a given
position. It filters the SNPs table for the given reference and
alternative alleles, and determines the genotype based on the allele
counts and percentages.

## Usage

``` r
sc_genotype(
  snps_tb,
  ref,
  alt,
  seqname,
  pos,
  alt_min_count = 2,
  alt_min_pct = 0.1,
  ref_min_count = 1,
  ref_min_pct = 1
)
```

## Arguments

- snps_tb:

  tibble: the SNPs table, output from `sc_mutations`.

- ref:

  character(1): the reference allele.

- alt:

  character(1): the alternative allele.

- seqname:

  character(1): the chromosome name of the position.

- pos:

  integer(1): the position of the mutation, 1-based.

- alt_min_count:

  integer(1): minimum UMI count of the alternative allele to call it
  "alt".

- alt_min_pct:

  numeric(1): minimum percentage of the alternative allele to call it
  "alt".

- ref_min_count:

  integer(1): minimum UMI count of the reference allele to call it
  "ref".

- ref_min_pct:

  numeric(1): minimum percentage of the reference allele to call it
  "ref".

## Value

A tibble with columns: barcode, allele_count_ref, pct_ref,
allele_count_alt, pct_alt, genotype.

## Examples

``` r
# get the SNPs table from sc_mutations
example(sc_mutations)
#> 
#> sc_mtt> ppl <- example_pipeline("SingleCellPipeline")
#> Writing configuration parameters to:  /tmp/RtmpmJ8vO7/file80d055968485/config_file_32976.json 
#> Warning: You have set to use oarfish quantification without gene quantification. Oarfish currently does not collapse UMIs, and gene quantification performs UMI collapsing. You may want to set do_gene_quantification to TRUE for more accurate results.
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: FALSE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> 
#> sc_mtt> ppl <- run_step(ppl, "barcode_demultiplex")
#> ── Running step: barcode_demultiplex @ Fri Oct 31 06:52:34 2025 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> FLEXIPLEX 0.96.2
#> Setting max barcode edit distance to 2
#> Setting max flanking sequence edit distance to 8
#> Setting read IDs to be  replaced
#> Setting number of threads to 8
#> Search pattern: 
#> primer: CTACACGACGCTCTTCCGATCT
#> BC: NNNNNNNNNNNNNNNN
#> UMI: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Setting known barcodes from /tmp/RtmpmJ8vO7/file80d055968485/bc_allow.tsv
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
#> 
#> sc_mtt> ppl <- run_step(ppl, "genome_alignment")
#> ── Running step: genome_alignment @ Fri Oct 31 06:52:34 2025 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpmJ8vO7/file80d055968485/matched_reads.fastq.gz -> /tmp/RtmpmJ8vO7/file80d055968485/align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> 
#> sc_mtt> snps_tb <- sc_mutations(
#> sc_mtt+   bam_path = ppl@genome_bam,
#> sc_mtt+   seqnames = c("chr14", "chr14"),
#> sc_mtt+   positions = c(1260, 2714), # positions of interest
#> sc_mtt+   indel = FALSE
#> sc_mtt+ )
#> 06:52:34 Got 1 bam file, parallelizing over each position ...
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%
#> 
#> 06:52:35 Merging results ...
#> 
#> sc_mtt> head(snps_tb)
#> # A tibble: 6 × 7
#>   allele barcode          allele_count cell_total_reads   pct   pos seqname
#>   <chr>  <chr>                   <dbl>            <dbl> <dbl> <dbl> <chr>  
#> 1 A      AACCATGAGTCGTTTG            0                2 0      1260 chr14  
#> 2 A      AACTCTTGTCACCTAA            0                1 0      1260 chr14  
#> 3 A      AACTTTCCACAGACTT            0                1 0      1260 chr14  
#> 4 A      AAGCCGCGTGTGAATA            0                4 0      1260 chr14  
#> 5 A      AAGGAGCGTGCTGTAT            1                3 0.333  1260 chr14  
#> 6 A      AATCGGTTCAGGTTCA            0                1 0      1260 chr14  
#> 
#> sc_mtt> snps_tb |>
#> sc_mtt+   dplyr::filter(pos == 1260) |>
#> sc_mtt+   dplyr::group_by(allele) |>
#> sc_mtt+   dplyr::summarise(count = sum(allele_count)) # should be identical to samtools pileup
#> # A tibble: 5 × 2
#>   allele count
#>   <chr>  <dbl>
#> 1 -         56
#> 2 A        103
#> 3 C          4
#> 4 G        169
#> 5 T          6
genotype_tb <- snps_tb |>
  sc_genotype(
    ref = "G", alt = "A", seqname = "chr14", pos = 1260,
    alt_min_count = 2, alt_min_pct = 0.1,
    ref_min_count = 1, ref_min_pct = 1
  )
dplyr::count(genotype_tb, genotype)
#> # A tibble: 3 × 2
#>   genotype     n
#>   <chr>    <int>
#> 1 alt         21
#> 2 ref         36
#> 3 NA          61
head(genotype_tb)
#> # A tibble: 6 × 6
#>   barcode          allele_count_alt allele_count_ref pct_alt pct_ref genotype
#>   <chr>                       <dbl>            <dbl>   <dbl>   <dbl> <chr>   
#> 1 AAGGAGCGTGCTGTAT                1                1   0.333   0.333 NA      
#> 2 ACCGTAACAAGCGTAG                2                3   0.4     0.6   alt     
#> 3 ACGATGTGTCTCCACT                2                0   1       0     alt     
#> 4 ACGCCAGAGAAGGCCT                1                0   1       0     NA      
#> 5 ACGCCGACAGATCCAT                1                0   1       0     NA      
#> 6 ACGCCGACATCACGTA                1                0   1       0     NA      
```
