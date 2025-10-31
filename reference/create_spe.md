# Create a SpatialExperiment object

This function creates a SpatialExperiment object from a
SingleCellExperiment object and a spatial barcode file.

## Usage

``` r
create_spe(
  sce,
  spatial_barcode_file,
  mannual_align_json,
  image,
  tissue_positions_file
)
```

## Arguments

- sce:

  The SingleCellExperiment object obtained from running the
  [`sc_long_pipeline`](https://mritchielab.github.io/FLAMES/reference/sc_long_pipeline.md)
  function.

- spatial_barcode_file:

  The path to the spatial barcode file, e.g.
  `"spaceranger-2.1.1/lib/python/cellranger/barcodes/visium-v2_coordinates.txt"`.

- mannual_align_json:

  The path to the mannual alignment json file.

- image:

  'DataFrame' containing the image data. See
  [`?SpatialExperiment::readImgData`](https://rdrr.io/pkg/SpatialExperiment/man/readImgData.html)
  and
  [`?SpatialExperiment::SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html).

- tissue_positions_file:

  The path to Visium positions file, e.g.
  `"spaceranger-2.1.1/lib/python/cellranger/barcodes/visium-v2_tissue_positions_list.csv"`.

## Value

A SpatialExperiment object.
