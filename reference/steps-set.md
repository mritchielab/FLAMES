# Set steps to perform in the pipeline

Set steps to perform in the pipeline

## Usage

``` r
steps(pipeline) <- value

# S4 method for class 'FLAMES.Pipeline'
steps(pipeline) <- value
```

## Arguments

- pipeline:

  An object of class \`FLAMES.Pipeline\`

- value:

  A named logical vector containing all possible steps for the pipeline.
  The names of the vector are the step names, and the values are logical
  indicating whether the step is configured to be performed.

## Value

An pipeline of class \`FLAMES.Pipeline\` with the updated steps.

## Examples

``` r
ppl <- example_pipeline()
#> Writing configuration parameters to:  /tmp/RtmpmJ8vO7/file80d03c791bb5/config_file_32976.json 
#> Warning: You have set to use oarfish quantification without gene quantification. Oarfish currently does not collapse UMIs, and gene quantification performs UMI collapsing. You may want to set do_gene_quantification to TRUE for more accurate results.
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: FALSE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
steps(ppl) <- c(
  barcode_demultiplex = TRUE,
  genome_alignment = TRUE,
  gene_quantification = TRUE,
  isoform_identification = FALSE,
  read_realignment = FALSE,
  transcript_quantification = TRUE
)
ppl
#> → A FLAMES.SingleCellPipeline outputting to /tmp/RtmpmJ8vO7/file80d03c791bb5
#> 
#> ── Inputs 
#> ✔ fastq: ...ibrary/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> ✔ annotation: /__w/_temp/Library/FLAMES/extdata/rps24.gtf.gz
#> ✔ genome_fa: /tmp/RtmpmJ8vO7/file80d03c791bb5/rps24.fa
#> ✔ barcodes_file: /tmp/RtmpmJ8vO7/file80d03c791bb5/bc_allow.tsv
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
#> ℹ transcript_quantification (pending)
# or partially change a step:
steps(ppl)["read_realignment"] <- TRUE
ppl
#> → A FLAMES.SingleCellPipeline outputting to /tmp/RtmpmJ8vO7/file80d03c791bb5
#> 
#> ── Inputs 
#> ✔ fastq: ...ibrary/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> ✔ annotation: /__w/_temp/Library/FLAMES/extdata/rps24.gtf.gz
#> ✔ genome_fa: /tmp/RtmpmJ8vO7/file80d03c791bb5/rps24.fa
#> ✔ barcodes_file: /tmp/RtmpmJ8vO7/file80d03c791bb5/bc_allow.tsv
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
#> ℹ read_realignment (pending)
#> ℹ transcript_quantification (pending)
```
