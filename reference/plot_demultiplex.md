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
#> ℹ Writing configuration to: /tmp/Rtmp4nGYdi/filebc447a72298d/config_file_48196.json
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> ── Running step: barcode_demultiplex @ Fri Jun 26 08:24:54 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/Rtmp4nGYdi/filebc447a72298d/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> CB:Z: tag field: CB
#> Processing file: /tmp/Rtmp4nGYdi/filebc447a72298d/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/Rtmp4nGYdi/filebc447a72298d/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/Rtmp4nGYdi/filebc447a72298d/fastq/sample3.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 993
#> Number of reads where at least one barcode was found: 932
#> Number of chimera reads: 2
#> All done!
#> Reads    Barcodes
#> 28   1
#> 27   1
#> 23   1
#> 22   1
#> 21   3
#> 20   2
#> 19   2
#> 18   1
#> 17   1
#> 15   1
#> 14   2
#> 13   3
#> 12   3
#> 11   6
#> 10   5
#> 9    6
#> 8    5
#> 7    5
#> 6    9
#> 5    14
#> 4    7
#> 3    36
#> 2    20
#> 1    2
#> Loading known barcodes from /tmp/Rtmp4nGYdi/filebc447a72298d/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> CB:Z: tag field: CB
#> Processing file: /tmp/Rtmp4nGYdi/filebc447a72298d/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 300
#> Number of reads where at least one barcode was found: 281
#> Number of chimera reads: 0
#> All done!
#> Reads    Barcodes
#> 9    1
#> 8    1
#> 7    5
#> 6    5
#> 5    2
#> 4    7
#> 3    18
#> 2    24
#> 1    62
#> Loading known barcodes from /tmp/Rtmp4nGYdi/filebc447a72298d/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> CB:Z: tag field: CB
#> Processing file: /tmp/Rtmp4nGYdi/filebc447a72298d/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 300
#> Number of reads where at least one barcode was found: 283
#> Number of chimera reads: 1
#> All done!
#> Reads    Barcodes
#> 9    2
#> 8    1
#> 7    2
#> 6    5
#> 5    3
#> 4    11
#> 3    17
#> 2    26
#> 1    53
#> Loading known barcodes from /tmp/Rtmp4nGYdi/filebc447a72298d/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> CB:Z: tag field: CB
#> Processing file: /tmp/Rtmp4nGYdi/filebc447a72298d/fastq/sample3.fq.gz
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
