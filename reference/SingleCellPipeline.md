# Pipeline for Single Cell Data

Semi-supervised isofrom detection and annotation for long read data.
This variant is meant for single sample scRNA-seq data. Specific
parameters can be configured in the config file (see
[`create_config`](https://mritchielab.github.io/FLAMES/reference/create_config.md)),
input files are specified via arguments.

## Usage

``` r
SingleCellPipeline(
  config_file,
  outdir,
  fastq,
  annotation,
  genome_fa,
  genome_mmi,
  minimap2,
  samtools,
  barcodes_file,
  expect_cell_number,
  controllers
)
```

## Arguments

- config_file:

  Path to the JSON configuration file. See
  [`create_config`](https://mritchielab.github.io/FLAMES/reference/create_config.md)
  for creating one.

- outdir:

  Path to the output directory. If it does not exist, it will be
  created.

- fastq:

  Path to the FASTQ file or a directory containing FASTQ files. Each
  file will be processed as an individual sample.

- annotation:

  The file path to the annotation file in GFF3 / GTF format.

- genome_fa:

  The file path to the reference genome in FASTA format.

- genome_mmi:

  (optional) The file path to minimap2's index reference genome.

- minimap2:

  (optional) The path to the minimap2 binary. If not provided, FLAMES
  will use a copy from bioconda via `basilisk`.

- samtools:

  (optional) The path to the samtools binary. If not provided, FLAMES
  will use a copy from bioconda via `basilisk`.

- barcodes_file:

  The file with expected cell barcodes, with each barcode on a new line.

- expect_cell_number:

  The expected number of cells in the sample. This is used if
  `barcodes_file` is not provided. See `BLAZE` for more details.

- controllers:

  (optional, **experimental**) A `crew_class_controller` object for
  running certain steps

## Value

A `FLAMES.SingleCellPipeline` object. The pipeline can be run using
`run_FLAMES(pipeline)`. The results can be accessed with
`experiment(pipeline)`. The pipeline also writes a number of files into
the given `outdir` directory, for example:

- matched_reads.fastq.gz:

  \- demultiplexed reads (barcode/UMI in the read header)

- align2genome.bam:

  \- sorted BAM file with reads aligned to the genome

- gene_count.mtx, gene_count_features.tsv, gene_count_barcodes.tsv:

  \- gene count matrix (Matrix Market format)

- isoform_annotated.gtf:

  \- updated annotation with novel isoforms (`.gff3` when not using
  bambu)

- transcript_assembly.fa:

  \- transcript sequences from the isoforms

- realign2transcript.bam:

  \- sorted realigned BAM file using transcript_assembly.fa as reference

- experiment.rds:

  \- the serialised SingleCellExperiment returned by
  `experiment(pipeline)`

See the *Expected output files* section of the FLAMES vignette
([`vignette("FLAMES_vignette")`](https://mritchielab.github.io/FLAMES/articles/FLAMES_vignette.md))
for the complete, per-step list (including the Oarfish and multi-sample
variants).

## Details

By default the pipeline starts with demultiplexing the input fastq data.
If the cell barcodes are known apriori (e.g. via coupled short-read
sequencing), the `barcodes_file` argument can be used to specify a file
containing the cell barcodes, and a modified Rcpp version of `flexiplex`
will be used; otherwise, `expect_cell_number` need to be provided, and
`BLAZE` will be used to generate the cell barcodes. The pipeline then
aligns the reads to the genome using `minimap2`. The alignment is then
used for isoform detection (either using `FLAMES` or `bambu`, can be
configured). The reads are then realigned to the detected isoforms.
Finally, a transcript count matrix is generated (either using `FLAMES`'s
simplistic counting or `oarfish`'s Expectation Maximization algorithm,
can be configured). The results can be accssed with
`experiment(pipeline)`. If the pipeline errored out / new steps were
configured, it can be resumed by calling `resume_FLAMES(pipeline)`

## See also

[`create_config`](https://mritchielab.github.io/FLAMES/reference/create_config.md)
for creating a configuration file,
[`BulkPipeline`](https://mritchielab.github.io/FLAMES/reference/BulkPipeline.md)
for bulk long data,
[`MultiSampleSCPipeline`](https://mritchielab.github.io/FLAMES/reference/MultiSampleSCPipeline.md)
for multi sample single cell pipelines.

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
ppl <- SingleCellPipeline(
  config_file = create_config(
    outdir,
    pipeline_parameters.demultiplexer = "flexiplex",
    pipeline_parameters.do_gene_quantification = FALSE
  ),
  outdir = outdir,
  fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
  annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  genome_fa = genome_fa,
  barcodes_file = bc_allow
)
#> ℹ Writing configuration to: /tmp/RtmpgehRKJ/filebbc02864eb80/config_file_48064.json
#> Warning: You have set to use oarfish quantification without gene quantification. Oarfish currently does not collapse UMIs, and gene quantification performs UMI collapsing. You may want to set do_gene_quantification to TRUE for more accurate results.
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: FALSE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
ppl <- run_FLAMES(ppl)
#> FLAMES version 2.7.1 (unknown source)
#> ── Running step: barcode_demultiplex @ Thu Aug 13 10:09:22 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/RtmpgehRKJ/filebbc02864eb80/bc_allow.tsv
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
#> ── Running step: genome_alignment @ Thu Aug 13 10:09:23 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpgehRKJ/filebbc02864eb80/matched_reads.fastq.gz -> /tmp/RtmpgehRKJ/filebbc02864eb80/align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: isoform_identification @ Thu Aug 13 10:09:23 2026 ─────────────
#> ── Running step: read_realignment @ Thu Aug 13 10:09:23 2026 ───────────────────
#> Checking for fastq file(s) /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpgehRKJ/filebbc02864eb80/matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpgehRKJ/filebbc02864eb80/matched_reads_dedup.fastq.gz
#>  files not found
#> Warning: Oarfish does not support UMI deduplication, you should deduplicate reads before running Oarfish
#> Realigning sample /tmp/RtmpgehRKJ/filebbc02864eb80/matched_reads.fastq.gz -> /tmp/RtmpgehRKJ/filebbc02864eb80/realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> ── Running step: transcript_quantification @ Thu Aug 13 10:09:23 2026 ──────────
experiment(ppl)
#> class: SingleCellExperiment 
#> dim: 10 137 
#> metadata(0):
#> assays(1): counts
#> rownames(10): ENSMUSG00000025290.17_19_5159_1
#>   ENSMUSG00000025290.17_19_5159_2 ... ENSMUST00000169826.2
#>   ENSMUST00000225023.1
#> rowData names(6): transcript_id source ... rank gene_id
#> colnames(137): AACCATGAGTCGTTTG AACTCTTGTCACCTAA ... TTGTAGGTCAGTGTTG
#>   TTTATGCAGACTAGAT
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```
