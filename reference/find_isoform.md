# Isoform identification

Long-read isoform identification with FLAMES or bambu.

## Usage

``` r
find_isoform(annotation, genome_fa, genome_bam, outdir, config)
```

## Arguments

- annotation:

  Path to annotation file. If configured to use bambu, the annotation
  must be provided as GTF file.

- genome_fa:

  The file path to genome fasta file.

- genome_bam:

  File path to BAM alignment file. Multiple files could be provided.

- outdir:

  The path to directory to store all output files.

- config:

  Parsed FLAMES configurations.

## Value

The updated annotation and the transcriptome assembly will be saved in
the output folder as `isoform_annotated.gff3` (GTF if bambu is selected)
and `transcript_assembly.fa` respectively.
