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
#> Writing configuration parameters to:  /tmp/Rtmpbzssfl/filea7097ada4e27/config_file_42761.json 
#> Warning: You have set to use oarfish quantification without gene quantification. Oarfish currently does not collapse UMIs, and gene quantification performs UMI collapsing. You may want to set do_gene_quantification to TRUE for more accurate results.
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: FALSE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> → A FLAMES.SingleCellPipeline outputting to /tmp/Rtmpbzssfl/filea7097ada4e27
#> 
#> ── Inputs 
#> ✔ fastq: ...ibrary/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> ✔ annotation: /__w/_temp/Library/FLAMES/extdata/rps24.gtf.gz
#> ✔ genome_fa: /tmp/Rtmpbzssfl/filea7097ada4e27/rps24.fa
#> ✔ barcodes_file: /tmp/Rtmpbzssfl/filea7097ada4e27/bc_allow.tsv
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
