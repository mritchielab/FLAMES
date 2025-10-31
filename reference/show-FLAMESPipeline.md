# Show method for FLAMES.Pipeline

Displays the pipeline in a pretty format

## Usage

``` r
# S4 method for class 'FLAMES.Pipeline'
show(object)

# S4 method for class 'FLAMES.SingleCellPipeline'
show(object)

# S4 method for class 'FLAMES.MultiSampleSCPipeline'
show(object)
```

## Arguments

- object:

  An object of class \`FLAMES.Pipeline\`

## Value

None. Displays output to the console.

## Examples

``` r
ppl <- example_pipeline()
#> Writing configuration parameters to:  /tmp/RtmpM0tt18/filea7097bfb3de4/config_file_42761.json 
#> Warning: You have set to use oarfish quantification without gene quantification. Oarfish currently does not collapse UMIs, and gene quantification performs UMI collapsing. You may want to set do_gene_quantification to TRUE for more accurate results.
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: FALSE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
show(ppl)
#> → A FLAMES.SingleCellPipeline outputting to /tmp/RtmpM0tt18/filea7097bfb3de4
#> 
#> ── Inputs 
#> ✔ fastq: ...ibrary/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> ✔ annotation: /__w/_temp/Library/FLAMES/extdata/rps24.gtf.gz
#> ✔ genome_fa: /tmp/RtmpM0tt18/filea7097bfb3de4/rps24.fa
#> ✔ barcodes_file: /tmp/RtmpM0tt18/filea7097bfb3de4/bc_allow.tsv
#> 
#> ── Outputs 
#> ℹ demultiplexed_fastq: matched_reads.fastq.gz
#> ℹ deduped_fastq: matched_reads_dedup.fastq.gz
#> ℹ genome_bam: align2genome.bam
#> ℹ transcriptome_assembly: transcript_assembly.fa
#> ℹ transcriptome_bam: realign2transcript.bam
#> 
#> ── Pipeline Steps 
#> ℹ barcode_demultiplex (pending)
#> ℹ genome_alignment (pending)
#> ℹ isoform_identification (pending)
#> ℹ read_realignment (pending)
#> ℹ transcript_quantification (pending)
```
