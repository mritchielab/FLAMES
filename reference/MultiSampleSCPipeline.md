# Pipeline for multi-sample long-read scRNA-seq data

Semi-supervised isofrom detection and annotation for long read data.
This variant is meant for multi-sample scRNA-seq data. Specific
parameters can be configured in the config file (see
[`create_config`](https://mritchielab.github.io/FLAMES/reference/create_config.md)),
input files are specified via arguments.

## Usage

``` r
MultiSampleSCPipeline(
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

  A named vector of fastq file (or folder) paths. Each element of the
  vector will be treated as a sample. The names of the vector will be
  used as the sample names. If not named, the sample names will be
  generated from the file names.

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

A `FLAMES.MultiSampleSCPipeline` object. The pipeline can be run using
the
[`run_FLAMES`](https://mritchielab.github.io/FLAMES/reference/run_FLAMES.md)
function. The resulting list of SingleCellExperiment objects can be
accessed using the `experiment` method.

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

[`SingleCellPipeline`](https://mritchielab.github.io/FLAMES/reference/SingleCellPipeline.md)
for single-sample long data and more details on the pipeline output,
[`create_config`](https://mritchielab.github.io/FLAMES/reference/create_config.md)
for creating a configuration file,
[`BulkPipeline`](https://mritchielab.github.io/FLAMES/reference/BulkPipeline.md)
for bulk long data.

## Examples

``` r
reads <- ShortRead::readFastq(
  system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES")
)
outdir <- tempfile()
dir.create(outdir)
dir.create(file.path(outdir, "fastq"))
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
ShortRead::writeFastq(reads[1:100],
  file.path(outdir, "fastq/sample1.fq.gz"), mode = "w", full = FALSE)
reads <- reads[-(1:100)]
ShortRead::writeFastq(reads[1:100],
  file.path(outdir, "fastq/sample2.fq.gz"), mode = "w", full = FALSE)
reads <- reads[-(1:100)]
ShortRead::writeFastq(reads,
  file.path(outdir, "fastq/sample3.fq.gz"), mode = "w", full = FALSE)
ppl <- MultiSampleSCPipeline(
  config_file = create_config(
    outdir,
    pipeline_parameters.demultiplexer = "flexiplex",
    pipeline_parameters.threads = 1,
    alignment_parameters.no_flank = TRUE
  ),
  outdir = outdir,
  fastq = c("sampleA" = file.path(outdir, "fastq"),
    "sample1" = file.path(outdir, "fastq", "sample1.fq.gz"),
    "sample2" = file.path(outdir, "fastq", "sample2.fq.gz"),
    "sample3" = file.path(outdir, "fastq", "sample3.fq.gz")),
  annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  genome_fa = genome_fa,
  barcodes_file = rep(bc_allow, 4)
)
#> Writing configuration parameters to:  /tmp/RtmpjRxi1E/file95d25647a0b/config_file_38354.json 
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
ppl <- run_FLAMES(ppl)
#> ── Running step: barcode_demultiplex @ Thu Feb 19 01:55:15 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/RtmpjRxi1E/file95d25647a0b/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpjRxi1E/file95d25647a0b/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/RtmpjRxi1E/file95d25647a0b/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/RtmpjRxi1E/file95d25647a0b/fastq/sample3.fq.gz
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
#> Loading known barcodes from /tmp/RtmpjRxi1E/file95d25647a0b/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpjRxi1E/file95d25647a0b/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 100
#> Number of reads where at least one barcode was found: 92
#> Number of chimera reads: 1
#> All done!
#> Reads    Barcodes
#> 4    1
#> 3    9
#> 2    9
#> 1    44
#> Loading known barcodes from /tmp/RtmpjRxi1E/file95d25647a0b/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpjRxi1E/file95d25647a0b/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 100
#> Number of reads where at least one barcode was found: 95
#> Number of chimera reads: 0
#> All done!
#> Reads    Barcodes
#> 4    2
#> 3    3
#> 2    16
#> 1    47
#> Loading known barcodes from /tmp/RtmpjRxi1E/file95d25647a0b/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpjRxi1E/file95d25647a0b/fastq/sample3.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 193
#> Number of reads where at least one barcode was found: 181
#> Number of chimera reads: 0
#> All done!
#> Reads    Barcodes
#> 7    1
#> 6    1
#> 5    1
#> 4    7
#> 3    10
#> 2    27
#> 1    53
#> ── Running step: genome_alignment @ Thu Feb 19 01:55:15 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpjRxi1E/file95d25647a0b/sampleA_matched_reads.fastq.gz -> /tmp/RtmpjRxi1E/file95d25647a0b/sampleA_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample /tmp/RtmpjRxi1E/file95d25647a0b/sample1_matched_reads.fastq.gz -> /tmp/RtmpjRxi1E/file95d25647a0b/sample1_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample /tmp/RtmpjRxi1E/file95d25647a0b/sample2_matched_reads.fastq.gz -> /tmp/RtmpjRxi1E/file95d25647a0b/sample2_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> Aligning sample /tmp/RtmpjRxi1E/file95d25647a0b/sample3_matched_reads.fastq.gz -> /tmp/RtmpjRxi1E/file95d25647a0b/sample3_align2genome.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 1 threads...
#> Indexing bam files
#> ── Running step: gene_quantification @ Thu Feb 19 01:55:16 2026 ────────────────
#> 01:55:16 AM Thu Feb 19 2026 quantify genes 
#> Using BAM(s): /tmp/RtmpjRxi1E/file95d25647a0b/sampleA_align2genome.bam,
#> /tmp/RtmpjRxi1E/file95d25647a0b/sample1_align2genome.bam,
#> /tmp/RtmpjRxi1E/file95d25647a0b/sample2_align2genome.bam, and
#> /tmp/RtmpjRxi1E/file95d25647a0b/sample3_align2genome.bam
#> ── Running step: isoform_identification @ Thu Feb 19 01:55:17 2026 ─────────────
#> ── Running step: read_realignment @ Thu Feb 19 01:55:18 2026 ───────────────────
#> Checking for fastq file(s) /tmp/RtmpjRxi1E/file95d25647a0b/fastq, /tmp/RtmpjRxi1E/file95d25647a0b/fastq/sample1.fq.gz, /tmp/RtmpjRxi1E/file95d25647a0b/fastq/sample2.fq.gz, /tmp/RtmpjRxi1E/file95d25647a0b/fastq/sample3.fq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpjRxi1E/file95d25647a0b/sampleA_matched_reads.fastq.gz, /tmp/RtmpjRxi1E/file95d25647a0b/sample1_matched_reads.fastq.gz, /tmp/RtmpjRxi1E/file95d25647a0b/sample2_matched_reads.fastq.gz, /tmp/RtmpjRxi1E/file95d25647a0b/sample3_matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpjRxi1E/file95d25647a0b/sampleA_matched_reads_dedup.fastq.gz, /tmp/RtmpjRxi1E/file95d25647a0b/sample1_matched_reads_dedup.fastq.gz, /tmp/RtmpjRxi1E/file95d25647a0b/sample2_matched_reads_dedup.fastq.gz, /tmp/RtmpjRxi1E/file95d25647a0b/sample3_matched_reads_dedup.fastq.gz
#>  files found
#> Realigning sample /tmp/RtmpjRxi1E/file95d25647a0b/sampleA_matched_reads_dedup.fastq.gz -> /tmp/RtmpjRxi1E/file95d25647a0b/sampleA_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 1 with CB threads...
#> Realigning sample /tmp/RtmpjRxi1E/file95d25647a0b/sample1_matched_reads_dedup.fastq.gz -> /tmp/RtmpjRxi1E/file95d25647a0b/sample1_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 1 with CB threads...
#> Realigning sample /tmp/RtmpjRxi1E/file95d25647a0b/sample2_matched_reads_dedup.fastq.gz -> /tmp/RtmpjRxi1E/file95d25647a0b/sample2_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 1 with CB threads...
#> Realigning sample /tmp/RtmpjRxi1E/file95d25647a0b/sample3_matched_reads_dedup.fastq.gz -> /tmp/RtmpjRxi1E/file95d25647a0b/sample3_realign2transcript.bam
#> Warning: samtools not found, using Rsamtools instead, this could be slower and might fail for large BAM files.
#> Sorting BAM files by 1 with CB threads...
#> ── Running step: transcript_quantification @ Thu Feb 19 01:55:19 2026 ──────────
experiment(ppl)
#> $sampleA
#> class: SingleCellExperiment 
#> dim: 14 137 
#> metadata(0):
#> assays(1): counts
#> rownames(14): ENSMUSG00000025290.17_19_5159_1
#>   ENSMUSG00000025290.17_19_5159_10 ... ENSMUST00000169826.2
#>   ENSMUST00000225023.1
#> rowData names(6): transcript_id source ... rank gene_id
#> colnames(137): CB:AACCATGAGTCGTTTG CB:AACTCTTGTCACCTAA ...
#>   CB:TTGTAGGTCTTACCTA CB:TTTATGCAGACTAGAT
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> 
#> $sample1
#> class: SingleCellExperiment 
#> dim: 14 62 
#> metadata(0):
#> assays(1): counts
#> rownames(14): ENSMUSG00000025290.17_19_5159_1
#>   ENSMUSG00000025290.17_19_5159_10 ... ENSMUST00000169826.2
#>   ENSMUST00000225023.1
#> rowData names(6): transcript_id source ... rank gene_id
#> colnames(62): CB:AAGCCGCGTGTGAATA CB:AAGGAGCGTGCTGTAT ...
#>   CB:TTGTAGGTCAGTGTTG CB:TTTATGCAGACTAGAT
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> 
#> $sample2
#> class: SingleCellExperiment 
#> dim: 14 67 
#> metadata(0):
#> assays(1): counts
#> rownames(14): ENSMUSG00000025290.17_19_5159_1
#>   ENSMUSG00000025290.17_19_5159_10 ... ENSMUST00000169826.2
#>   ENSMUST00000225023.1
#> rowData names(6): transcript_id source ... rank gene_id
#> colnames(67): CB:AACCATGAGTCGTTTG CB:AAGCCGCGTGTGAATA ...
#>   CB:TTGTAGGTCTTACCTA CB:TTTATGCAGACTAGAT
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> 
#> $sample3
#> class: SingleCellExperiment 
#> dim: 14 99 
#> metadata(0):
#> assays(1): counts
#> rownames(14): ENSMUSG00000025290.17_19_5159_1
#>   ENSMUSG00000025290.17_19_5159_10 ... ENSMUST00000169826.2
#>   ENSMUST00000225023.1
#> rowData names(6): transcript_id source ... rank gene_id
#> colnames(99): CB:AACTCTTGTCACCTAA CB:AACTTTCCACAGACTT ...
#>   CB:TTGTAGGGTGGGTATG CB:TTGTAGGTCAGTGTTG
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> 
```
