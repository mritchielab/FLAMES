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
#> Writing configuration parameters to:  /tmp/RtmpmJ8vO7/file80d0586f6d37/config_file_32976.json 
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
#> ── Running step: barcode_demultiplex @ Fri Oct 31 06:51:23 2025 ────────────────
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
#> Setting known barcodes from /tmp/RtmpmJ8vO7/file80d0586f6d37/bc_allow.tsv
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
#> ── Running step: genome_alignment @ Fri Oct 31 06:51:23 2025 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpmJ8vO7/file80d0586f6d37/matched_reads.fastq.gz -> /tmp/RtmpmJ8vO7/file80d0586f6d37/align2genome.bam
#> Your fastq file appears to have tags, but you did not provide the -y option to minimap2 to include the tags in the output.
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: gene_quantification @ Fri Oct 31 06:51:24 2025 ────────────────
#> 06:51:24 AM Fri Oct 31 2025 quantify genes 
#> Using BAM(s): /tmp/RtmpmJ8vO7/file80d0586f6d37/align2genome.bam
#> ── Running step: isoform_identification @ Fri Oct 31 06:51:25 2025 ─────────────
#> ── Running step: read_realignment @ Fri Oct 31 06:51:25 2025 ───────────────────
#> Checking for fastq file(s) /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpmJ8vO7/file80d0586f6d37/matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpmJ8vO7/file80d0586f6d37/matched_reads_dedup.fastq.gz
#>  files found
#> Realigning sample /tmp/RtmpmJ8vO7/file80d0586f6d37/matched_reads_dedup.fastq.gz -> /tmp/RtmpmJ8vO7/file80d0586f6d37/realign2transcript.bam
#> Your fastq file appears to have tags, but you did not provide the -y option to minimap2 to include the tags in the output.
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: transcript_quantification @ Fri Oct 31 06:51:26 2025 ──────────
#> 06:51:26 AM Fri Oct 31 2025 quantify transcripts 
#> Warning: Annotation in GFF format may cause errors. Please consider using GTF formats.
#> Found realignment file(s):   realign2transcript.bam
#> Pipeline saved to /tmp/RtmpmJ8vO7/file80d0586f6d37/pipeline.rds
sce_2 <- create_sce_from_dir(outdir, annotation)
```
