# Plot spatial pie chart of isoforms

This function plots a spatial pie chart for given features.

## Usage

``` r
plot_spatial_isoform(spe, isoforms, assay_type = "counts", color_palette, ...)
```

## Arguments

- spe:

  The SpatialExperiment object.

- isoforms:

  The isoforms to plot.

- assay_type:

  The assay that contains the given features. E.g. 'counts',
  'logcounts'.

- color_palette:

  Named vector of colors for each isoform.

- ...:

  Additional arguments to pass to
  [`plot_spatial_pie`](https://mritchielab.github.io/FLAMES/reference/plot_spatial_pie.md),
  including `opacity`, `grayscale`, `pie_scale`.

## Value

A ggplot object.
