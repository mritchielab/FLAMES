# Variant count for single-cell data

Count the number of reads supporting each variants at the given
positions for each cell.

## Usage

``` r
sc_mutations(bam_path, seqnames, positions, indel = FALSE, threads = 1)
```

## Arguments

- bam_path:

  character(1) or character(n): path to the bam file(s) aligned to the
  reference genome (NOT the transcriptome! Unless the postions are also
  from the transcriptome).

- seqnames:

  character(n): chromosome names of the postions to count alleles.

- positions:

  integer(n): positions, 1-based, same length as seqnames. The positions
  to count alleles.

- indel:

  logical(1): whether to count indels (TRUE) or SNPs (FALSE).

- threads:

  integer(1): number of threads to use. Maximum number of threads is the
  number of bam files \* number of positions.

## Value

A tibble with columns: allele, barcode, allele_count, cell_total_reads,
pct, pos, seqname.

## Examples

``` r
ppl <- example_pipeline("SingleCellPipeline")
#> Writing configuration parameters to:  /tmp/RtmpaNKtnp/file80d04d33c591/config_file_32976.json 
#> Warning: You have set to use oarfish quantification without gene quantification. Oarfish currently does not collapse UMIs, and gene quantification performs UMI collapsing. You may want to set do_gene_quantification to TRUE for more accurate results.
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: FALSE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
ppl <- run_step(ppl, "barcode_demultiplex")
#> ── Running step: barcode_demultiplex @ Fri Oct 31 06:44:33 2025 ────────────────
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
#> Setting known barcodes from /tmp/RtmpaNKtnp/file80d04d33c591/bc_allow.tsv
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
ppl <- run_step(ppl, "genome_alignment")
#> ── Running step: genome_alignment @ Fri Oct 31 06:44:34 2025 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpaNKtnp/file80d04d33c591/matched_reads.fastq.gz -> /tmp/RtmpaNKtnp/file80d04d33c591/align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
snps_tb <- sc_mutations(
  bam_path = ppl@genome_bam,
  seqnames = c("chr14", "chr14"),
  positions = c(1260, 2714), # positions of interest
  indel = FALSE
)
#> 06:44:34 Got 1 bam file, parallelizing over each position ...
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%
#> 
#> 06:44:35 Merging results ...
head(snps_tb)
#> # A tibble: 6 × 7
#>   allele barcode          allele_count cell_total_reads   pct   pos seqname
#>   <chr>  <chr>                   <dbl>            <dbl> <dbl> <dbl> <chr>  
#> 1 A      AACCATGAGTCGTTTG            0                2 0      1260 chr14  
#> 2 A      AACTCTTGTCACCTAA            0                1 0      1260 chr14  
#> 3 A      AACTTTCCACAGACTT            0                1 0      1260 chr14  
#> 4 A      AAGCCGCGTGTGAATA            0                4 0      1260 chr14  
#> 5 A      AAGGAGCGTGCTGTAT            1                3 0.333  1260 chr14  
#> 6 A      AATCGGTTCAGGTTCA            0                1 0      1260 chr14  
snps_tb |>
  dplyr::filter(pos == 1260) |>
  dplyr::group_by(allele) |>
  dplyr::summarise(count = sum(allele_count)) # should be identical to samtools pileup
#> # A tibble: 5 × 2
#>   allele count
#>   <chr>  <dbl>
#> 1 -         56
#> 2 A        103
#> 3 C          4
#> 4 G        169
#> 5 T          6
```
