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
#> ℹ Writing configuration to: /tmp/Rtmp4nGYdi/filebc442c7b4a4e/config_file_48196.json
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
show(ppl)
#> → A FLAMES.SingleCellPipeline outputting to /tmp/Rtmp4nGYdi/filebc442c7b4a4e
#> 
#> ── Inputs 
#> ✔ fastq: ...ibrary/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> ✔ annotation: /__w/_temp/Library/FLAMES/extdata/rps24.gtf.gz
#> ✔ genome_fa: /tmp/Rtmp4nGYdi/filebc442c7b4a4e/rps24.fa
#> ✔ barcodes_file: /tmp/Rtmp4nGYdi/filebc442c7b4a4e/bc_allow.tsv
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
#> ℹ gene_quantification (pending)
#> ℹ isoform_identification (pending)
#> ℹ read_realignment (pending)
#> ℹ transcript_quantification (pending)
```
