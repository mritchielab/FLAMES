# Plot spatial pie chart

This function plots a spatial pie chart for given features.

## Usage

``` r
plot_spatial_pie(
  spe,
  features,
  assay_type = "counts",
  color_palette,
  opacity = 50,
  grayscale = TRUE,
  pie_scale = 0.8
)
```

## Arguments

- spe:

  The SpatialExperiment object.

- features:

  The features to plot.

- assay_type:

  The assay that contains the given features.

- color_palette:

  Named vector of colors for each feature.

- opacity:

  The opacity of the background tissue image.

- grayscale:

  Whether to convert the background image to grayscale.

- pie_scale:

  The size of the pie charts.

## Value

A ggplot object.
