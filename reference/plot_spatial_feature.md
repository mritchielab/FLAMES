# Plot feature on spatial image

This function plots a spatial point plot for given feature

## Usage

``` r
plot_spatial_feature(
  spe,
  feature,
  opacity = 50,
  grayscale = TRUE,
  size = 1,
  assay_type = "counts",
  color = "red",
  ...
)
```

## Arguments

- spe:

  The SpatialExperiment object.

- feature:

  The feature to plot. Could be either a feature name or index present
  in the assay or a numeric vector of length nrow(spe).

- opacity:

  The opacity of the background tissue image.

- grayscale:

  Whether to convert the background image to grayscale.

- size:

  The size of the points.

- assay_type:

  The assay that contains the given features. E.g. 'counts',
  'logcounts'.

- color:

  The maximum color for the feature. Minimum color is transparent.

- ...:

  Additional arguments to pass to
  [`geom_point`](https://ggplot2.tidyverse.org/reference/geom_point.html).

## Value

A ggplot object.
