# Create `SingleCellExperiment` object from `FLAMES` output folder

Create `SingleCellExperiment` object from `FLAMES` output folder

## Usage

``` r
create_sce_from_dir(outdir, annotation, quantification = "FLAMES")
```

## Arguments

- outdir:

  The folder containing `FLAMES` output files

- annotation:

  the annotation file that was used to produce the output files

- quantification:

  (Optional) the quantification method used to generate the output files
  (either "FLAMES" or "Oarfish".). If not specified, the function will
  attempt to determine the quantification method.

## Value

a list of `SingleCellExperiment` objects if multiple transcript matrices
were found in the output folder, or a `SingleCellExperiment` object if
only one were found

## Examples

``` r
outdir <- tempfile()
dir.create(outdir)
bc_allow <- file.path(outdir, "bc_allow.tsv")
genome_fa <- file.path(outdir, "rps24.fa")
R.utils::gunzip(
  filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
  destname = bc_allow, remove = FALSE
)
R.utils::gunzip(
  filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
  destname = genome_fa, remove = FALSE
)
annotation <- system.file("extdata", "rps24.gtf.gz", package = "FLAMES")

sce <- sc_long_pipeline(
  genome_fa = genome_fa,
  fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
  annotation = annotation,
  outdir = outdir,
  barcodes_file = bc_allow,
  config_file = create_config(
    outdir,
    pipeline_parameters.demultiplexer = "flexiplex",
    oarfish_quantification = FALSE
  )
)
#> ℹ Writing configuration to: /tmp/RtmpnC89xy/filebc1553cd34fa/config_file_48149.json
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> ── Running step: barcode_demultiplex @ Thu Aug  6 04:48:29 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/RtmpnC89xy/filebc1553cd34fa/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 8
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> CB:Z: tag field: CB
#> Processing file: /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
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
#> ── Running step: genome_alignment @ Thu Aug  6 04:48:29 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpnC89xy/filebc1553cd34fa/matched_reads.fastq.gz -> /tmp/RtmpnC89xy/filebc1553cd34fa/align2genome.bam
#> Your fastq file appears to have tags, but you did not provide the -y option to minimap2 to include the tags in the output.
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: gene_quantification @ Thu Aug  6 04:48:29 2026 ────────────────
#> 04:48:29 AM Thu Aug 06 2026 quantify genes 
#> Using BAM(s): /tmp/RtmpnC89xy/filebc1553cd34fa/align2genome.bam
#> ── Running step: isoform_identification @ Thu Aug  6 04:48:30 2026 ─────────────
#> ── Running step: read_realignment @ Thu Aug  6 04:48:31 2026 ───────────────────
#> Checking for fastq file(s) /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpnC89xy/filebc1553cd34fa/matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpnC89xy/filebc1553cd34fa/matched_reads_dedup.fastq.gz
#>  files found
#> Realigning sample /tmp/RtmpnC89xy/filebc1553cd34fa/matched_reads_dedup.fastq.gz -> /tmp/RtmpnC89xy/filebc1553cd34fa/realign2transcript.bam
#> Your fastq file appears to have tags, but you did not provide the -y option to minimap2 to include the tags in the output.
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: transcript_quantification @ Thu Aug  6 04:48:31 2026 ──────────
#> 04:48:31 AM Thu Aug 06 2026 quantify transcripts 
#> Warning: Annotation in GFF format may cause errors. Please consider using GTF formats.
#> Found realignment file(s):   realign2transcript.bam
#> Pipeline saved to /tmp/RtmpnC89xy/filebc1553cd34fa/pipeline.rds
sce_2 <- create_sce_from_dir(outdir, annotation)
```
