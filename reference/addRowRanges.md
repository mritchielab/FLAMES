# Add rowRanges by rownames to `SummarizedExperiment` object Assumes rownames are transcript_ids Assumes transcript_id is present in the annotation file

Add rowRanges by rownames to `SummarizedExperiment` object Assumes
rownames are transcript_ids Assumes transcript_id is present in the
annotation file

## Usage

``` r
addRowRanges(sce, annotation, outdir)
```

## Value

a `SummarizedExperiment` object with rowRanges added
