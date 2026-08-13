# Get pipeline configurations

This function returns the configuration of the pipeline.

## Usage

``` r
config(pipeline)

# S4 method for class 'FLAMES.Pipeline'
config(pipeline)
```

## Arguments

- pipeline:

  An object of class \`FLAMES.Pipeline\`.

## Value

A list containing the configuration of the pipeline.

## Examples

``` r
pipeline <- example_pipeline(type = "BulkPipeline")
#> ℹ Writing configuration to: /tmp/RtmpgehRKJ/filebbc02c05b13d/config_file_48064.json
#> Configured steps: 
#>  genome_alignment: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
config(pipeline)
#> $pipeline_parameters
#> $pipeline_parameters$seed
#> [1] 2022
#> 
#> $pipeline_parameters$threads
#> [1] 1
#> 
#> $pipeline_parameters$demultiplexer
#> [1] "BLAZE"
#> 
#> $pipeline_parameters$do_barcode_demultiplex
#> [1] TRUE
#> 
#> $pipeline_parameters$do_gene_quantification
#> [1] TRUE
#> 
#> $pipeline_parameters$do_genome_alignment
#> [1] TRUE
#> 
#> $pipeline_parameters$do_isoform_identification
#> [1] TRUE
#> 
#> $pipeline_parameters$bambu_isoform_identification
#> [1] FALSE
#> 
#> $pipeline_parameters$multithread_isoform_identification
#> [1] FALSE
#> 
#> $pipeline_parameters$do_read_realignment
#> [1] TRUE
#> 
#> $pipeline_parameters$do_transcript_quantification
#> [1] TRUE
#> 
#> $pipeline_parameters$oarfish_quantification
#> [1] TRUE
#> 
#> 
#> $additional_arguments
#> $additional_arguments$oarfish
#> list()
#> 
#> $additional_arguments$BLAZE
#> [1] "--10x-kit-version" "3v3"              
#> 
#> 
#> $barcode_parameters
#> $barcode_parameters$max_bc_editdistance
#> [1] 2
#> 
#> $barcode_parameters$max_flank_editdistance
#> [1] 8
#> 
#> $barcode_parameters$segments
#>      type                pattern   name bc_list_name group buffer_size
#> 1   FIXED CTACACGACGCTCTTCCGATCT primer         <NA>  <NA>          NA
#> 2 MATCHED       NNNNNNNNNNNNNNNN     CB           CB                 5
#> 3  RANDOM           NNNNNNNNNNNN     UB         <NA>  <NA>          NA
#> 4   FIXED              TTTTTTTTT  polyT         <NA>  <NA>          NA
#>   max_edit_distance
#> 1                NA
#> 2                 2
#> 3                NA
#> 4                NA
#> 
#> $barcode_parameters$barcode_groups
#> list()
#> 
#> $barcode_parameters$strand
#> [1] "-"
#> 
#> $barcode_parameters$TSO_seq
#> [1] "AAGCAGTGGTATCAACGCAGAGTACATGGG"
#> 
#> $barcode_parameters$TSO_prime
#> [1] 5
#> 
#> $barcode_parameters$cutadapt_minimum_length
#> [1] 10
#> 
#> $barcode_parameters$full_length_only
#> [1] FALSE
#> 
#> 
#> $isoform_parameters
#> $isoform_parameters$generate_raw_isoform
#> [1] FALSE
#> 
#> $isoform_parameters$max_dist
#> [1] 10
#> 
#> $isoform_parameters$max_ts_dist
#> [1] 100
#> 
#> $isoform_parameters$max_splice_match_dist
#> [1] 10
#> 
#> $isoform_parameters$min_fl_exon_len
#> [1] 40
#> 
#> $isoform_parameters$max_site_per_splice
#> [1] 3
#> 
#> $isoform_parameters$min_sup_cnt
#> [1] 5
#> 
#> $isoform_parameters$min_cnt_pct
#> [1] 0.001
#> 
#> $isoform_parameters$min_sup_pct
#> [1] 0.2
#> 
#> $isoform_parameters$bambu_ndr
#> [1] 0.5
#> 
#> $isoform_parameters$bambu_verbose
#> [1] FALSE
#> 
#> $isoform_parameters$bambu_trust_reference
#> [1] TRUE
#> 
#> $isoform_parameters$strand_specific
#> [1] 0
#> 
#> $isoform_parameters$remove_incomp_reads
#> [1] 4
#> 
#> $isoform_parameters$downsample_ratio
#> [1] 1
#> 
#> 
#> $alignment_parameters
#> $alignment_parameters$use_junctions
#> [1] TRUE
#> 
#> $alignment_parameters$no_flank
#> [1] TRUE
#> 
#> 
#> $realign_parameters
#> $realign_parameters$use_annotation
#> [1] TRUE
#> 
#> 
#> $transcript_counting
#> $transcript_counting$min_tr_coverage
#> [1] 0.4
#> 
#> $transcript_counting$min_read_coverage
#> [1] 0.4
#> 
#> 
```
