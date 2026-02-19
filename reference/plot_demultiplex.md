# Plot Cell Barcode demultiplex statistics

produce a barplot of cell barcode demultiplex statistics

## Usage

``` r
plot_demultiplex(pipeline)

# S4 method for class 'FLAMES.SingleCellPipeline'
plot_demultiplex(pipeline)
```

## Arguments

- pipeline:

  A `FLAMES.SingleCellPipeline` object

## Value

a list of ggplot objects:

- reads_count_plot: stacked barplot of: demultiplexed reads

- knee_plot: knee plot of UMI counts before TSO trimming

- flank_editdistance_plot: flanking sequence (adaptor) edit-distance
  plot

- barcode_editdistance_plot: barcode edit-distance plot

- cutadapt_plot: if TSO trimming is performed, number of reads kept by
  cutadapt

## Examples

``` r
pipeline <- example_pipeline("MultiSampleSCPipeline") |>
  run_step("barcode_demultiplex")
#> Writing configuration parameters to:  /tmp/RtmpjRxi1E/file95d26f6ab3bf/config_file_38354.json 
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> ── Running step: barcode_demultiplex @ Thu Feb 19 01:56:05 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/RtmpjRxi1E/file95d26f6ab3bf/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpjRxi1E/file95d26f6ab3bf/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/RtmpjRxi1E/file95d26f6ab3bf/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/RtmpjRxi1E/file95d26f6ab3bf/fastq/sample3.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 993
#> Number of reads where at least one barcode was found: 928
#> Number of chimera reads: 3
#> All done!
#> Reads    Barcodes
#> 28   1
#> 26   1
#> 24   1
#> 21   1
#> 20   2
#> 19   3
#> 18   1
#> 17   3
#> 16   2
#> 14   2
#> 13   1
#> 12   6
#> 11   4
#> 10   4
#> 9    7
#> 8    5
#> 7    4
#> 6    13
#> 5    14
#> 4    3
#> 3    39
#> 2    15
#> 1    5
#> Loading known barcodes from /tmp/RtmpjRxi1E/file95d26f6ab3bf/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpjRxi1E/file95d26f6ab3bf/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 300
#> Number of reads where at least one barcode was found: 279
#> Number of chimera reads: 1
#> All done!
#> Reads    Barcodes
#> 9    2
#> 7    2
#> 6    2
#> 5    8
#> 4    9
#> 3    16
#> 2    31
#> 1    52
#> Loading known barcodes from /tmp/RtmpjRxi1E/file95d26f6ab3bf/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpjRxi1E/file95d26f6ab3bf/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 300
#> Number of reads where at least one barcode was found: 281
#> Number of chimera reads: 1
#> All done!
#> Reads    Barcodes
#> 9    1
#> 8    2
#> 7    3
#> 6    3
#> 5    5
#> 4    11
#> 3    14
#> 2    22
#> 1    64
#> Loading known barcodes from /tmp/RtmpjRxi1E/file95d26f6ab3bf/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpjRxi1E/file95d26f6ab3bf/fastq/sample3.fq.gz
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
plot_demultiplex(pipeline)
#> $reads_count_plot

#> 
#> $knee_plot
#> `geom_smooth()` using formula = 'y ~ x'

#> 
#> $flank_editdistance_plot

#> 
#> $cutadapt_plot

#> 
```
