# FLAMES 2.3.1

## FLAMES

The overhauled FLAMES 2.3.1 pipeline provides convenient pipelines for
performing single-cell and bulk long-read analysis of mutations and
isoforms. The pipeline is designed to take various type of experiments,
e.g. with or without known cell barcodes and custome cell barcode
designs, to reduce the need of constantly re-inventing the wheel for
every new sequencing protocol.

![(#fig:workflow) FLAMES pipeline](FLAMESpipeline-01.png)

(#fig:workflow) FLAMES pipeline

### Creating a pipeline

To start your long-read RNA-seq analysis, simply create a pipeline via
either
[`BulkPipeline()`](https://mritchielab.github.io/FLAMES/reference/BulkPipeline.md),
[`SingleCellPipeline()`](https://mritchielab.github.io/FLAMES/reference/SingleCellPipeline.md)
or
[`MultiSampleSCPipeline()`](https://mritchielab.github.io/FLAMES/reference/MultiSampleSCPipeline.md).
Let’s try
[`SingleCellPipeline()`](https://mritchielab.github.io/FLAMES/reference/SingleCellPipeline.md)
first:

``` r

outdir <- tempfile()
dir.create(outdir)
# some example data
# known cell barcodes, e.g. from coupled short-read sequencing
bc_allow <- file.path(outdir, "bc_allow.tsv")
R.utils::gunzip(
  filename = system.file("extdata", "bc_allow.tsv.gz", package = "FLAMES"),
  destname = bc_allow, remove = FALSE
)
# reference genome
genome_fa <- file.path(outdir, "rps24.fa")
R.utils::gunzip(
  filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
  destname = genome_fa, remove = FALSE
)

pipeline <- SingleCellPipeline(
  # use the default configs
  config_file = create_config(
    outdir,
    pipeline_parameters.demultiplexer = "flexiplex"
  ),
  outdir = outdir,
  # the input fastq file
  fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
  # reference annotation file
  annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  genome_fa = genome_fa,
  barcodes_file = bc_allow
)
#> Writing configuration parameters to:  /tmp/RtmpkJkeWg/filea0cd3863f954/config_file_41165.json
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
pipeline
#> → A FLAMES.SingleCellPipeline outputting to /tmp/RtmpkJkeWg/filea0cd3863f954
#> 
#> ── Inputs
#> ✔ fastq: ...ibrary/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> ✔ annotation: /__w/_temp/Library/FLAMES/extdata/rps24.gtf.gz
#> ✔ genome_fa: /tmp/RtmpkJkeWg/filea0cd3863f954/rps24.fa
#> ✔ barcodes_file: /tmp/RtmpkJkeWg/filea0cd3863f954/bc_allow.tsv
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

### Running the pipeline

To execute the pipeline, simply call `run_FLAMES(pipeline)`. This will
run through all the steps in the pipeline, returning a updated pipeline
object:

``` r

pipeline <- run_FLAMES(pipeline)
#> ── Running step: barcode_demultiplex @ Thu Feb 19 01:57:25 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/RtmpkJkeWg/filea0cd3863f954/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 8
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
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
#> ── Running step: genome_alignment @ Thu Feb 19 01:57:25 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpkJkeWg/filea0cd3863f954/matched_reads.fastq.gz -> /tmp/RtmpkJkeWg/filea0cd3863f954/align2genome.bam
#> Warning in minimap2_align(fq_in = fastqs[i], fa_file = genome, config =
#> pipeline@config, : samtools not found, using Rsamtools instead, this could be
#> slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: gene_quantification @ Thu Feb 19 01:57:26 2026 ────────────────
#> 01:57:26 AM Thu Feb 19 2026 quantify genes
#> Using BAM(s): /tmp/RtmpkJkeWg/filea0cd3863f954/align2genome.bam
#> Assigning reads to genes...
#> Writing the gene count matrix ...
#> Plotting the saturation curve ...
#> Generating deduplicated fastq file ...
#> ── Running step: isoform_identification @ Thu Feb 19 01:57:26 2026 ─────────────
#> #### Read gene annotations
#>  Removed similar transcripts in gene annotation: Counter()
#> #### find isoforms
#> chr14
#> ── Running step: read_realignment @ Thu Feb 19 01:57:27 2026 ───────────────────
#> Checking for fastq file(s) /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpkJkeWg/filea0cd3863f954/matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpkJkeWg/filea0cd3863f954/matched_reads_dedup.fastq.gz
#>  files found
#> Realigning sample /tmp/RtmpkJkeWg/filea0cd3863f954/matched_reads_dedup.fastq.gz -> /tmp/RtmpkJkeWg/filea0cd3863f954/realign2transcript.bam
#> Warning in minimap2_align(fq_in = fastqs[i], fa_file =
#> pipeline@transcriptome_assembly, : samtools not found, using Rsamtools instead,
#> this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> ── Running step: transcript_quantification @ Thu Feb 19 01:57:27 2026 ──────────
pipeline
#> ✔ A FLAMES.SingleCellPipeline outputting to /tmp/RtmpkJkeWg/filea0cd3863f954
#> 
#> ── Inputs
#> ✔ fastq: ...ibrary/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> ✔ annotation: /__w/_temp/Library/FLAMES/extdata/rps24.gtf.gz
#> ✔ genome_fa: /tmp/RtmpkJkeWg/filea0cd3863f954/rps24.fa
#> ✔ barcodes_file: /tmp/RtmpkJkeWg/filea0cd3863f954/bc_allow.tsv
#> 
#> ── Outputs
#> ✔ demultiplexed_fastq: matched_reads.fastq.gz [219.7 KB]
#> ✔ deduped_fastq: matched_reads_dedup.fastq.gz [206.1 KB]
#> ✔ genome_bam: align2genome.bam [273.8 KB]
#> ✔ novel_isoform_annotation: isoform_annotated.gff3 [7.4 KB]
#> ✔ transcriptome_assembly: transcript_assembly.fa [8.4 KB]
#> ✔ transcriptome_bam: realign2transcript.bam [395.4 KB]
#> 
#> ── Pipeline Steps
#> ✔ barcode_demultiplex (completed in 0.40 sec)
#> ✔ genome_alignment (completed in 0.38 sec)
#> ✔ gene_quantification (completed in 0.74 sec)
#> ✔ isoform_identification (completed in 0.45 sec)
#> ✔ read_realignment (completed in 0.28 sec)
#> ✔ transcript_quantification (completed in 0.49 sec)
```

If you run into any error,
[`run_FLAMES()`](https://mritchielab.github.io/FLAMES/reference/run_FLAMES.md)
will stop and return the pipeline object with the error message. After
resolving the error, you can run `resume_FLAMES(pipeline)` to continue
the pipeline from the last step. There is also
`run_step(pipeline, step_name)` to run a specific step in the pipeline.
Let’s show this by deliberately causing an error via deleting the input
files:

``` r

# set up a new pipeline
outdir2 <- tempfile()
pipeline2 <- SingleCellPipeline(
  config_file = create_config(
    outdir,
    pipeline_parameters.demultiplexer = "flexiplex"
  ),
  outdir = outdir2,
  fastq = system.file("extdata", "fastq", "musc_rps24.fastq.gz", package = "FLAMES"),
  annotation = system.file("extdata", "rps24.gtf.gz", package = "FLAMES"),
  genome_fa = genome_fa,
  barcodes_file = bc_allow
)
#> Output directory does not exists: one is being created
#> [1] "/tmp/RtmpkJkeWg/filea0cd723a6cf1"
#> Writing configuration parameters to:  /tmp/RtmpkJkeWg/filea0cd3863f954/config_file_41165.json
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead

# delete the reference genome
unlink(genome_fa)
pipeline2 <- run_FLAMES(pipeline2)
#> ── Running step: barcode_demultiplex @ Thu Feb 19 01:57:28 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/RtmpkJkeWg/filea0cd3863f954/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 8
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
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
#> ── Running step: genome_alignment @ Thu Feb 19 01:57:28 2026 ───────────────────
#> Creating junction bed file from GFF3 annotation.
#> Aligning sample /tmp/RtmpkJkeWg/filea0cd723a6cf1/matched_reads.fastq.gz -> /tmp/RtmpkJkeWg/filea0cd723a6cf1/align2genome.bam
#> Warning in minimap2_align(fq_in = fastqs[i], fa_file = genome, config =
#> pipeline@config, : samtools not found, using Rsamtools instead, this could be
#> slower and might fail for large BAM files.
#> Warning in check_status_code(minimap2_status, cmd, "Minimap2"):
#> "'/__w/_temp/Library/FLAMES/bin/minimap2' '-ax' 'splice' '-k14'
#> '--secondary=no' '-t' '8' '--seed' '2022' '-y' '--junc-bed'
#> '/tmp/RtmpkJkeWg/filea0cd723a6cf1/reference.bed' '--junc-bonus' '1'
#> '/tmp/RtmpkJkeWg/filea0cd3863f954/rps24.fa'
#> '/tmp/RtmpkJkeWg/filea0cd723a6cf1/matched_reads.fastq.gz' >
#> '/tmp/RtmpkJkeWg/filea0cd723a6cf1/filea0cd49995913.sam'" exited with status
#> code 1.
#> Warning in value[[3L]](cond): Error in step genome_alignment: argument "no" is
#> missing, with no default, pipeline stopped.
pipeline2
#> ! A FLAMES.SingleCellPipeline outputting to /tmp/RtmpkJkeWg/filea0cd723a6cf1
#> 
#> ── Inputs
#> ✔ fastq: ...ibrary/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> ✔ annotation: /__w/_temp/Library/FLAMES/extdata/rps24.gtf.gz
#> ! genome_fa: /tmp/RtmpkJkeWg/filea0cd3863f954/rps24.fa [missing]
#> ✔ barcodes_file: /tmp/RtmpkJkeWg/filea0cd3863f954/bc_allow.tsv
#> 
#> ── Outputs
#> ✔ demultiplexed_fastq: matched_reads.fastq.gz [219.7 KB]
#> ℹ deduped_fastq: matched_reads_dedup.fastq.gz
#> ℹ genome_bam: align2genome.bam
#> ℹ transcriptome_assembly: transcript_assembly.fa
#> ℹ transcriptome_bam: realign2transcript.bam
#> 
#> ── Pipeline Steps
#> ✔ barcode_demultiplex (completed in 0.17 sec)
#> ✖ genome_alignment (failed: Error in ifelse(status_code == 137, "This is likely due to running out of memory."): argument "no" is missing, with no default
#> )
#> ℹ gene_quantification (pending)
#> ℹ isoform_identification (pending)
#> ℹ read_realignment (pending)
#> ℹ transcript_quantification (pending)
```

Let’s then fix the error by re-creating the reference genome file and
resume the pipeline:

``` r

R.utils::gunzip(
  filename = system.file("extdata", "rps24.fa.gz", package = "FLAMES"),
  destname = genome_fa, remove = FALSE
)
pipeline2 <- resume_FLAMES(pipeline2)
#> Resuming pipeline from step: genome_alignment
#> ── Running step: genome_alignment @ Thu Feb 19 01:57:29 2026 ───────────────────
#> Aligning sample /tmp/RtmpkJkeWg/filea0cd723a6cf1/matched_reads.fastq.gz -> /tmp/RtmpkJkeWg/filea0cd723a6cf1/align2genome.bam
#> Warning in minimap2_align(fq_in = fastqs[i], fa_file = genome, config =
#> pipeline@config, : samtools not found, using Rsamtools instead, this could be
#> slower and might fail for large BAM files.
#> Sorting BAM files by genome coordinates with 8 threads...
#> Indexing bam files
#> ── Running step: gene_quantification @ Thu Feb 19 01:57:29 2026 ────────────────
#> 01:57:29 AM Thu Feb 19 2026 quantify genes
#> Using BAM(s): /tmp/RtmpkJkeWg/filea0cd723a6cf1/align2genome.bam
#> Assigning reads to genes...
#> Writing the gene count matrix ...
#> Plotting the saturation curve ...
#> Generating deduplicated fastq file ...
#> ── Running step: isoform_identification @ Thu Feb 19 01:57:30 2026 ─────────────
#> #### Read gene annotations
#>  Removed similar transcripts in gene annotation: Counter()
#> #### find isoforms
#> chr14
#> ── Running step: read_realignment @ Thu Feb 19 01:57:30 2026 ───────────────────
#> Checking for fastq file(s) /__w/_temp/Library/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpkJkeWg/filea0cd723a6cf1/matched_reads.fastq.gz
#>  files found
#> Checking for fastq file(s) /tmp/RtmpkJkeWg/filea0cd723a6cf1/matched_reads_dedup.fastq.gz
#>  files found
#> Realigning sample /tmp/RtmpkJkeWg/filea0cd723a6cf1/matched_reads_dedup.fastq.gz -> /tmp/RtmpkJkeWg/filea0cd723a6cf1/realign2transcript.bam
#> Warning in minimap2_align(fq_in = fastqs[i], fa_file =
#> pipeline@transcriptome_assembly, : samtools not found, using Rsamtools instead,
#> this could be slower and might fail for large BAM files.
#> Sorting BAM files by 8 with CB threads...
#> ── Running step: transcript_quantification @ Thu Feb 19 01:57:30 2026 ──────────
pipeline2
#> ✔ A FLAMES.SingleCellPipeline outputting to /tmp/RtmpkJkeWg/filea0cd723a6cf1
#> 
#> ── Inputs
#> ✔ fastq: ...ibrary/FLAMES/extdata/fastq/musc_rps24.fastq.gz
#> ✔ annotation: /__w/_temp/Library/FLAMES/extdata/rps24.gtf.gz
#> ✔ genome_fa: /tmp/RtmpkJkeWg/filea0cd3863f954/rps24.fa
#> ✔ barcodes_file: /tmp/RtmpkJkeWg/filea0cd3863f954/bc_allow.tsv
#> 
#> ── Outputs
#> ✔ demultiplexed_fastq: matched_reads.fastq.gz [219.7 KB]
#> ✔ deduped_fastq: matched_reads_dedup.fastq.gz [206.1 KB]
#> ✔ genome_bam: align2genome.bam [273.8 KB]
#> ✔ novel_isoform_annotation: isoform_annotated.gff3 [7.4 KB]
#> ✔ transcriptome_assembly: transcript_assembly.fa [8.4 KB]
#> ✔ transcriptome_bam: realign2transcript.bam [395.4 KB]
#> 
#> ── Pipeline Steps
#> ✔ barcode_demultiplex (completed in 0.17 sec)
#> ✔ genome_alignment (completed in 0.26 sec)
#> ✔ gene_quantification (completed in 0.36 sec)
#> ✔ isoform_identification (completed in 0.26 sec)
#> ✔ read_realignment (completed in 0.27 sec)
#> ✔ transcript_quantification (completed in 0.42 sec)
```

After completing the pipeline, a `SingleCellExperiment` object is
created (or `SummarizedExperiment` for bulk pipeline and list of
`SingleCellExperiment` for multi-sample pipeline). You can access the
results via `experiment(pipeline)`:

``` r

experiment(pipeline)
#> class: SingleCellExperiment 
#> dim: 10 137 
#> metadata(0):
#> assays(1): counts
#> rownames(10): ENSMUSG00000025290.17_19_5159_1
#>   ENSMUSG00000025290.17_19_5159_2 ... ENSMUST00000169826.2
#>   ENSMUST00000225023.1
#> rowData names(6): transcript_id source ... rank gene_id
#> colnames(137): CB:AACCATGAGTCGTTTG CB:AACTCTTGTCACCTAA ...
#>   CB:TTGTAGGTCAGTGTTG CB:TTTATGCAGACTAGAT
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

#### HPC support

Individual steps can be submitted as HPC jobs via `crew` and
`crew.cluster` packages, simply supply a list of crew controllers (named
by step name) to the `controllers` argument of the pipeline
constructors. For example, we could run the alignment steps through
controllers, while keeping the rest in the main R session.

``` r

# example_pipeline provides an example pipeline for each of the three types
# of pipelines: BulkPipeline, SingleCellPipeline and MultiSampleSCPipeline
mspipeline <- example_pipeline("MultiSampleSCPipeline")
#> Writing configuration parameters to:  /tmp/RtmpkJkeWg/filea0cd18bfd67a/config_file_41165.json
#> Configured steps: 
#>  barcode_demultiplex: TRUE
#>  genome_alignment: TRUE
#>  gene_quantification: TRUE
#>  isoform_identification: TRUE
#>  read_realignment: TRUE
#>  transcript_quantification: TRUE
#> samtools not found, will use Rsamtools package instead
# Providing a single controller will run all steps in it:
controllers(mspipeline) <- crew::crew_controller_local()
# Setting controllers to an empty list will run all steps in the main R session:
controllers(mspipeline) <- list()
# Alternatively, we can run only the alignment steps in the crew controller:
controllers(mspipeline)[["genome_alignment"]] <- crew::crew_controller_local(workers = 4)
# Or `controllers(mspipeline) <- list(genome_alignment = crew::crew_cluster())`
# to remove controllers for all other steps.
# Replace `crew_controller_local()` with `crew.cluster::crew_controller_slurm()` or other
# crew controllers according to your HPC job scheduler.
```

[`run_FLAMES()`](https://mritchielab.github.io/FLAMES/reference/run_FLAMES.md)
will then submit the alignment step to the crew controller. By default,
[`run_step()`](https://mritchielab.github.io/FLAMES/reference/run_step.md)
will ignore the crew controllers and run the step in the main R session,
since it is easier to debug. You can use
`run_step(pipeline, step_name, disable_controller = FALSE)` to prevent
this behavior and run the step in crew controllers if available.

You can tailor the resources for each step by specifying different
arguments to the controllers. The alignment step typically benifits from
more cores (e.g. 64 cores and 20GB memory), whereas other steps might
not need as much cores but more memory hungry.

``` r

# An example helper function to create a Slurm controller with specific resources
create_slurm_controller <- function(
    cpus, memory_gb, workers = 10, seconds_idle = 10,
    script_lines = "module load R/flexiblas/4.5.0") {
  name <- sprintf("slurm_%dc%dg", cpus, memory_gb)
  crew.cluster::crew_controller_slurm(
    name = name,
    workers = workers,
    seconds_idle = seconds_idle,
    retry_tasks = FALSE,
    options_cluster = crew.cluster::crew_options_slurm(
      script_lines = script_lines,
      memory_gigabytes_required = memory_gb,
      cpus_per_task = cpus,
      log_output = file.path("logs", "crew_log_%A.txt"),
      log_error = file.path("logs", "crew_log_%A.txt")
    )
  )
}
controllers(mspipeline)[["genome_alignment"]] <-
  create_slurm_controller(cpus = 64, memory_gb = 20)
```

See also the [`crew.cluster`
website](https://wlandau.github.io/crew.cluster/reference/index.html)
for more information on the supported job schedulers.

### Visualizations

#### QC plots

Quality metrics are collected throughout the pipeline, and FLAMES
provide visiualization functions to plot the metrics. For the first
demultiplexing step, we can use the `plot_demultiplex` function to see
how well many reads are retained after demultiplexing:

``` r

# don't have to run the entire pipeline for this
# let's just run the demultiplexing step
mspipeline <- run_step(mspipeline, "barcode_demultiplex")
#> ── Running step: barcode_demultiplex @ Thu Feb 19 01:57:31 2026 ────────────────
#> Using flexiplex for barcode demultiplexing.
#> Loading known barcodes from /tmp/RtmpkJkeWg/filea0cd18bfd67a/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpkJkeWg/filea0cd18bfd67a/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/RtmpkJkeWg/filea0cd18bfd67a/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Processing file: /tmp/RtmpkJkeWg/filea0cd18bfd67a/fastq/sample3.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 993
#> Number of reads where at least one barcode was found: 929
#> Number of chimera reads: 2
#> All done!
#> Reads    Barcodes
#> 26   1
#> 25   1
#> 22   1
#> 21   1
#> 20   4
#> 19   1
#> 18   2
#> 17   3
#> 16   1
#> 15   3
#> 14   3
#> 12   2
#> 11   5
#> 10   5
#> 9    8
#> 8    4
#> 7    4
#> 6    9
#> 5    15
#> 4    6
#> 3    29
#> 2    26
#> 1    3
#> Loading known barcodes from /tmp/RtmpkJkeWg/filea0cd18bfd67a/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpkJkeWg/filea0cd18bfd67a/fastq/sample1.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 300
#> Number of reads where at least one barcode was found: 280
#> Number of chimera reads: 1
#> All done!
#> Reads    Barcodes
#> 9    1
#> 8    1
#> 7    1
#> 6    6
#> 5    10
#> 4    5
#> 3    14
#> 2    28
#> 1    56
#> Loading known barcodes from /tmp/RtmpkJkeWg/filea0cd18bfd67a/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpkJkeWg/filea0cd18bfd67a/fastq/sample2.fq.gz
#> Searching for barcodes...
#> Number of reads processed: 300
#> Number of reads where at least one barcode was found: 281
#> Number of chimera reads: 0
#> All done!
#> Reads    Barcodes
#> 8    1
#> 7    3
#> 6    5
#> 5    9
#> 4    11
#> 3    13
#> 2    20
#> 1    56
#> Loading known barcodes from /tmp/RtmpkJkeWg/filea0cd18bfd67a/bc_allow.tsv
#> Number of known barcodes: 143
#> FLEXIPLEX 1.02.6
#> Setting max flanking sequence edit distance to 8
#> Setting number of threads to 1
#> Search pattern:
#> primer: CTACACGACGCTCTTCCGATCT
#> CB: NNNNNNNNNNNNNNNN
#> UB: NNNNNNNNNNNN
#> polyT: TTTTTTTTT
#> Processing file: /tmp/RtmpkJkeWg/filea0cd18bfd67a/fastq/sample3.fq.gz
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
plot_demultiplex(mspipeline)
#> $reads_count_plot
```

![](FLAMES_vignette_files/figure-html/unnamed-chunk-8-1.png)

    #> 
    #> $knee_plot
    #> `geom_smooth()` using formula = 'y ~ x'

![](FLAMES_vignette_files/figure-html/unnamed-chunk-8-2.png)

    #> 
    #> $flank_editdistance_plot

![](FLAMES_vignette_files/figure-html/unnamed-chunk-8-3.png)

    #> 
    #> $cutadapt_plot

![](FLAMES_vignette_files/figure-html/unnamed-chunk-8-4.png)

#### Work in progress

More examples coming soon.

### FLAMES on Windows

Due to FLAMES requiring minimap2 and pysam, FLAMES is currently
unavaliable on Windows.

### Citation

Please cite the flames’s paper (Tian et al. 2020) if you use flames in
your research. As FLAMES used incorporates BLAZE (You et al. 2023),
flexiplex (Davidson et al. 2023) and minimap2 (Li 2018), samtools, bambu
(Chen et al. 2023). Please make sure to cite when using these tools.

## Session Info

    #> R Under development (unstable) (2026-02-14 r89420)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.3 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    #>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    #>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    #>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    #> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    #> 
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] FLAMES_2.5.2     BiocStyle_2.39.0
    #> 
    #> loaded via a namespace (and not attached):
    #>   [1] splines_4.6.0               later_1.4.6                
    #>   [3] BiocIO_1.21.0               bitops_1.0-9               
    #>   [5] filelock_1.0.3              tibble_3.3.1               
    #>   [7] R.oo_1.27.1                 polyclip_1.10-7            
    #>   [9] bambu_3.13.1                XML_3.99-0.22              
    #>  [11] lifecycle_1.0.5             pwalign_1.7.0              
    #>  [13] edgeR_4.9.2                 doParallel_1.0.17          
    #>  [15] vroom_1.7.0                 processx_3.8.6             
    #>  [17] lattice_0.22-9              MASS_7.3-65                
    #>  [19] magrittr_2.0.4              limma_3.67.0               
    #>  [21] sass_0.4.10                 rmarkdown_2.30             
    #>  [23] jquerylib_0.1.4             yaml_2.3.12                
    #>  [25] metapod_1.19.1              otel_0.2.0                 
    #>  [27] reticulate_1.45.0           cowplot_1.2.0              
    #>  [29] DBI_1.2.3                   RColorBrewer_1.1-3         
    #>  [31] abind_1.4-8                 ShortRead_1.69.2           
    #>  [33] GenomicRanges_1.63.1        purrr_1.2.1                
    #>  [35] R.utils_2.13.0              BiocGenerics_0.57.0        
    #>  [37] RCurl_1.98-1.17             yulab.utils_0.2.4          
    #>  [39] tweenr_2.0.3                rappdirs_0.3.4             
    #>  [41] circlize_0.4.17             IRanges_2.45.0             
    #>  [43] S4Vectors_0.49.0            ggrepel_0.9.6              
    #>  [45] irlba_2.3.7                 dqrng_0.4.1                
    #>  [47] pkgdown_2.2.0.9000          codetools_0.2-20           
    #>  [49] DelayedArray_0.37.0         scuttle_1.21.0             
    #>  [51] ggforce_0.5.0               tidyselect_1.2.1           
    #>  [53] shape_1.4.6.1               UCSC.utils_1.7.1           
    #>  [55] farver_2.1.2                ScaledMatrix_1.19.0        
    #>  [57] viridis_0.6.5               matrixStats_1.5.0          
    #>  [59] stats4_4.6.0                Seqinfo_1.1.0              
    #>  [61] GenomicAlignments_1.47.0    jsonlite_2.0.0             
    #>  [63] GetoptLong_1.1.0            BiocNeighbors_2.5.4        
    #>  [65] scater_1.39.2               iterators_1.0.14           
    #>  [67] systemfonts_1.3.1           foreach_1.5.2              
    #>  [69] tools_4.6.0                 ragg_1.5.0                 
    #>  [71] collections_0.3.11          Rcpp_1.1.1                 
    #>  [73] glue_1.8.0                  gridExtra_2.3              
    #>  [75] SparseArray_1.11.10         mgcv_1.9-4                 
    #>  [77] xfun_0.56                   MatrixGenerics_1.23.0      
    #>  [79] GenomeInfoDb_1.47.2         dplyr_1.2.0                
    #>  [81] withr_3.0.2                 BiocManager_1.30.27        
    #>  [83] fastmap_1.2.0               basilisk_1.23.0            
    #>  [85] bluster_1.21.0              latticeExtra_0.6-31        
    #>  [87] digest_0.6.39               rsvd_1.0.5                 
    #>  [89] R6_2.6.1                    textshaping_1.0.4          
    #>  [91] colorspace_2.1-2            jpeg_0.1-11                
    #>  [93] dichromat_2.0-0.1           RSQLite_2.4.6              
    #>  [95] cigarillo_1.1.0             R.methodsS3_1.8.2          
    #>  [97] tidyr_1.3.2                 generics_0.1.4             
    #>  [99] data.table_1.18.2.1         rtracklayer_1.71.3         
    #> [101] httr_1.4.8                  htmlwidgets_1.6.4          
    #> [103] S4Arrays_1.11.1             scatterpie_0.2.6           
    #> [105] pkgconfig_2.0.3             gtable_0.3.6               
    #> [107] blob_1.3.0                  ComplexHeatmap_2.27.1      
    #> [109] S7_0.2.1                    hwriter_1.3.2.1            
    #> [111] SingleCellExperiment_1.33.0 XVector_0.51.0             
    #> [113] htmltools_0.5.9             bookdown_0.46              
    #> [115] clue_0.3-67                 scales_1.4.0               
    #> [117] Biobase_2.71.0              png_0.1-8                  
    #> [119] nanonext_1.8.0              SpatialExperiment_1.21.0   
    #> [121] scran_1.39.0                ggfun_0.2.0                
    #> [123] knitr_1.51                  tzdb_0.5.0                 
    #> [125] rjson_0.2.23                nlme_3.1-168               
    #> [127] curl_7.0.0                  crew_1.3.0                 
    #> [129] cachem_1.1.0                GlobalOptions_0.1.3        
    #> [131] stringr_1.6.0               parallel_4.6.0             
    #> [133] vipor_0.4.7                 AnnotationDbi_1.73.0       
    #> [135] restfulr_0.0.16             desc_1.4.3                 
    #> [137] pillar_1.11.1               grid_4.6.0                 
    #> [139] vctrs_0.7.1                 promises_1.5.0             
    #> [141] BiocSingular_1.27.1         beachmat_2.27.2            
    #> [143] cluster_2.1.8.2             beeswarm_0.4.0             
    #> [145] evaluate_1.0.5              readr_2.1.6                
    #> [147] GenomicFeatures_1.63.1      magick_2.9.0               
    #> [149] locfit_1.5-9.12             cli_3.6.5                  
    #> [151] compiler_4.6.0              Rsamtools_2.27.0           
    #> [153] rlang_1.1.7                 crayon_1.5.3               
    #> [155] labeling_0.4.3              interp_1.1-6               
    #> [157] ps_1.9.1                    fs_1.6.6                   
    #> [159] ggbeeswarm_0.7.3            stringi_1.8.7              
    #> [161] viridisLite_0.4.3           deldir_2.0-4               
    #> [163] BiocParallel_1.45.0         Biostrings_2.79.4          
    #> [165] Matrix_1.7-4                dir.expiry_1.19.0          
    #> [167] BSgenome_1.79.1             hms_1.1.4                  
    #> [169] bit64_4.6.0-1               ggplot2_4.0.2              
    #> [171] statmod_1.5.1               KEGGREST_1.51.1            
    #> [173] SummarizedExperiment_1.41.1 mirai_2.6.0                
    #> [175] igraph_2.2.2                memoise_2.0.1              
    #> [177] bslib_0.10.0                bit_4.6.0                  
    #> [179] xgboost_3.2.0.1

## References

Chen, Ying, Andre Sim, Yuk Kei Wan, et al. 2023. “Context-Aware
Transcript Quantification from Long-Read RNA-Seq Data with Bambu.”
*Nature Methods*, 1–9.

Davidson, Nadia M, Noorul Amin, Ling Min Hao, et al. 2023. “Flexiplex: A
Versatile Demultiplexer and Search Tool for Omics Data.” *bioRxiv*,
2023–08.

Li, Heng. 2018. “Minimap2: pairwise alignment for nucleotide sequences.”
*Bioinformatics* 34 (18): 3094–100.
<https://doi.org/10.1093/bioinformatics/bty191>.

Tian, Luyi, Jafar S. Jabbari, Rachel Thijssen, et al. 2020.
“Comprehensive Characterization of Single Cell Full-Length Isoforms in
Human and Mouse with Long-Read Sequencing.” *bioRxiv*, ahead of print.
<https://doi.org/10.1101/2020.08.10.243543>.

You, Yupei, Yair DJ Prawer, De Paoli-Iseppi, et al. 2023.
“Identification of Cell Barcodes from Long-Read Single-Cell RNA-Seq with
BLAZE.” *Genome Biology* 24 (1): 1–23.
