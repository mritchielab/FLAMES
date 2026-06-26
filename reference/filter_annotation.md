# filter annotation for plotting coverages

Removes isoform annotations that could produce ambigious reads, such as
isoforms that only differ by the 5' / 3' end. This could be useful for
plotting average coverage plots.

## Usage

``` r
filter_annotation(annotation, keep = "tss_differ")
```

## Arguments

- annotation:

  path to the GTF annotation file, or the parsed GenomicRanges object
  with a valid `transcript_id` column, and each Range representing a
  transcript.

- keep:

  string, one of 'tss_differ' (only keep isoforms that all differ by the
  transcription start site position), 'tes_differ' (only keep those that
  differ by the transcription end site position), 'both' (only keep
  those that differ by both the start and end site), or
  'single_transcripts' (only keep genes that contain a single
  transcript).

## Value

GenomicRanges of the filtered isoforms

## Examples

``` r
filtered_annotation <- filter_annotation(
  system.file("extdata", "rps24.gtf.gz", package = 'FLAMES'), keep = 'tes_differ')
filtered_annotation
#> GRanges object with 6 ranges and 18 metadata columns:
#>       seqnames    ranges strand |   source       type     score     phase
#>          <Rle> <IRanges>  <Rle> | <factor>   <factor> <numeric> <integer>
#>   [1]    chr14   19-5159      + |   HAVANA transcript        NA      <NA>
#>   [2]    chr14   32-3389      + |   HAVANA transcript        NA      <NA>
#>   [3]    chr14   68-5124      + |   HAVANA transcript        NA      <NA>
#>   [4]    chr14   86-1118      + |   HAVANA transcript        NA      <NA>
#>   [5]    chr14  160-2761      + |   HAVANA transcript        NA      <NA>
#>   [6]    chr14  450-1290      + |   HAVANA transcript        NA      <NA>
#>                     gene_id        transcript_id      gene_type   gene_name
#>                 <character>          <character>    <character> <character>
#>   [1] ENSMUSG00000025290.17 ENSMUST00000225994.1 protein_coding       Rps24
#>   [2] ENSMUSG00000025290.17 ENSMUST00000225117.1 protein_coding       Rps24
#>   [3] ENSMUSG00000025290.17 ENSMUST00000224568.1 protein_coding       Rps24
#>   [4] ENSMUSG00000025290.17 ENSMUST00000224549.1 protein_coding       Rps24
#>   [5] ENSMUSG00000025290.17 ENSMUST00000224569.1 protein_coding       Rps24
#>   [6] ENSMUSG00000025290.17 ENSMUST00000224699.1 protein_coding       Rps24
#>            transcript_type transcript_name       level           protein_id
#>                <character>     <character> <character>          <character>
#>   [1]      retained_intron       Rps24-212           2                 <NA>
#>   [2] processed_transcript       Rps24-211           2                 <NA>
#>   [3]       protein_coding       Rps24-207           2 ENSMUSP00000153637.1
#>   [4] processed_transcript       Rps24-206           2                 <NA>
#>   [5] processed_transcript       Rps24-208           2                 <NA>
#>   [6] processed_transcript       Rps24-209           2                 <NA>
#>       transcript_support_level      mgi_id         tag      ccdsid
#>                    <character> <character> <character> <character>
#>   [1]                     <NA>   MGI:98147        <NA>        <NA>
#>   [2]                     <NA>   MGI:98147        <NA>        <NA>
#>   [3]                     <NA>   MGI:98147       basic        <NA>
#>   [4]                     <NA>   MGI:98147 mRNA_end_NF        <NA>
#>   [5]                     <NA>   MGI:98147        <NA>        <NA>
#>   [6]                     <NA>   MGI:98147 mRNA_end_NF        <NA>
#>                havana_gene    havana_transcript
#>                <character>          <character>
#>   [1] OTTMUSG00000068639.1 OTTMUST00000165910.1
#>   [2] OTTMUSG00000068639.1 OTTMUST00000165914.1
#>   [3] OTTMUSG00000068639.1 OTTMUST00000165915.1
#>   [4] OTTMUSG00000068639.1 OTTMUST00000165916.1
#>   [5] OTTMUSG00000068639.1 OTTMUST00000165917.1
#>   [6] OTTMUSG00000068639.1 OTTMUST00000165918.1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
