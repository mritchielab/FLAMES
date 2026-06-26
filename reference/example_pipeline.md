# Example pipelins

Provides example pipelines for bulk, single cell and multi-sample single
cell.

## Usage

``` r
example_pipeline(type = "SingleCellPipeline", outdir)
```

## Arguments

- type:

  The type of pipeline to create. Options are "SingleCellPipeline",
  "BulkPipeline", and "MultiSampleSCPipeline".

- outdir:

  (Optional) The output directory where the example pipeline will be
  created. If not provided, a temporary directory will be created.

## Value

A pipeline object of the specified type.

## See also

[`SingleCellPipeline`](https://mritchielab.github.io/FLAMES/reference/SingleCellPipeline.md)
for creating the single cell pipeline,
[`BulkPipeline`](https://mritchielab.github.io/FLAMES/reference/BulkPipeline.md)
for bulk long data,
[`MultiSampleSCPipeline`](https://mritchielab.github.io/FLAMES/reference/MultiSampleSCPipeline.md)
for multi sample single cell pipelines.

## Examples

``` r
example_pipeline("SingleCellPipeline")
#> ℹ Writing configuration to: /tmp/Rtmp4nGYdi/filebc44b24d89c/config_file_48196.json
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> → A FLAMES.SingleCellPipeline outputting to /tmp/Rtmp4nGYdi/filebc44b24d89c
#> 
#> ── Inputs 
#> ✔ fastq: ...ibrary/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> ✔ annotation: /__w/_temp/Library/FLAMES/extdata/rps24.gtf.gz
#> ✔ genome_fa: /tmp/Rtmp4nGYdi/filebc44b24d89c/rps24.fa
#> ✔ barcodes_file: /tmp/Rtmp4nGYdi/filebc44b24d89c/bc_allow.tsv
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
