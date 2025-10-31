# Gene quantification

Calculate the per gene UMI count matrix by parsing the genome alignment
file.

## Usage

``` r
quantify_gene(
  annotation,
  outdir,
  pipeline = "sc_single_sample",
  infq,
  in_bam,
  out_fastq,
  n_process,
  saturation_curve = TRUE,
  sample_names = NULL,
  random_seed = 2024
)
```

## Arguments

- annotation:

  The file path to the annotation file in GFF3 format

- outdir:

  The path to directory to store all output files.

- pipeline:

  The pipeline type as a character string, either `sc_single_sample`
  (single-cell, single-sample), `bulk` (bulk, single or multi-sample),
  or `sc_multi_sample` (single-cell, multiple samples)

- infq:

  The input FASTQ file.

- in_bam:

  The input BAM file(s) from the genome alignment step.

- out_fastq:

  The output FASTQ file(s) to store deduplicated reads.

- n_process:

  The number of processes to use for parallelization.

- saturation_curve:

  Logical, whether to generate a saturation curve figure.

- sample_names:

  A vector of sample names, default to the file names of input fastq
  files, or folder names if `fastqs` is a vector of folders.

- random_seed:

  The random seed for reproducibility.

## Value

The count matrix will be saved in the output folder as
`transcript_count.csv.gz`.

## Details

After the genome alignment step (`do_genome_align`), the alignment file
will be parsed to generate the per gene UMI count matrix. For each gene
in the annotation file, the number of reads overlapping with the gene’s
genomic coordinates will be assigned to that gene. If a read overlaps
multiple genes, it will be assigned to the gene with the highest number
of overlapping nucleotides. If exon coordinates are included in the
provided annotation, the decision will first consider the number of
nucleotides aligned to the exons of each gene. In cases of a tie, the
overlap with introns will be used as a tiebreaker. If there is still a
tie after considering both exons and introns, a random gene will be
selected from the tied candidates.

After the read-to-gene assignment, the per gene UMI count matrix will be
generated. Specifically, for each gene, the reads with similar mapping
coordinates of transcript termination sites (TTS, i.e. the end of the
the read with a polyT or polyA) will be grouped together. UMIs of reads
in the same group will be collapsed to generate the UMI counts for each
gene.

Finally, a new fastq file with deduplicated reads by keeping the longest
read in each UMI.
