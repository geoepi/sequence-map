<p align="center">
  <img src="images/sequence-map-sticker.svg" width="220" alt="sequence-map hex sticker">
</p>

# sequencemap

`sequencemap` is the R package used by the `sequence-map` repository.
It converts aligned viral sequences and sampling coordinates into PCA-first,
location-level summaries and optional spatial predictions of ordination-axis
scores across geographic space.

## What sequence-map does

`sequencemap` supports:

- FASTA alignment input and sequence-metadata validation.
- Alignment QC and one-hot encoding of canonical sequence variation.
- PCA-first ordination of sequence variation.
- Optional DAPC when meaningful group labels are supplied.
- Location-level score and reduced-space dispersion summaries.
- Projected prediction grids and INLA/SPDE spatial axis models.
- Raster GeoTIFF exports of predicted axis surfaces.
- Polygon/province summaries derived from rasterized prediction surfaces.

## What sequence-map does not do

The package does **not** estimate nucleotide diversity directly from predicted
surfaces. PCA/LD axis surfaces are predicted genetic-position scores. Polygon
summaries are reporting products derived from continuous model predictions, not
separate province-level models or direct nucleotide-diversity estimates.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("JMHumphreys/sequence-map")
```

Or, from a local clone of the repository:

```r
install.packages("devtools")
devtools::install(".")
```

For development work, use:

```r
devtools::load_all()
devtools::document()
testthat::test_dir("tests/testthat")
```

## INLA installation

`INLA` is listed under `Suggests` because installation is nonstandard and is
only required for mesh construction, SPDE axis-model fitting, and surface
prediction. Non-spatial preprocessing, PCA, DAPC, and location summaries can run
without INLA.

For the stable INLA release, use the official INLA repository:

```r
install.packages(
  "INLA",
  repos = c(
    getOption("repos"),
    INLA = "https://inla.r-inla-download.org/R/stable"
  ),
  dep = TRUE
)
```

To update an existing INLA installation:

```r
INLA::inla.upgrade()
```

Windows users may need to uninstall INLA before reinstalling. Some INLA
workflows may also require Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("graph", "Rgraphviz"), dep = TRUE)
```

The INLA project notes that most installation problems are caused by outdated R
or INLA versions, and recommends using a recent R version and the current INLA
repository release.

## Inputs

A typical workflow requires:

- A FASTA alignment with unique sequence identifiers.
- Sequence metadata with matching identifiers and numeric longitude/latitude.
- A boundary readable by `sf` for prediction-grid construction.
- A projected CRS in metric units for spatial modeling.

## Core workflow

```text
read_alignment()
-> validate_alignment_metadata_match()
-> alignment_to_variant_matrix()
-> run_sequence_pca()
-> summarize_location_scores()
-> summarize_location_diversity()
-> prepare_spatial_points()
-> make_prediction_grid()
-> make_inla_mesh()
-> fit_axis_models()
-> predict_axis_surfaces()
-> axis_surface_to_raster()
-> aggregate_axis_rasters_to_polygons()
-> calculate_polygon_axis_metrics()
```

DAPC is opt-in because it requires a biologically justified grouping decision or
explicitly exploratory inferred groups.

## Minimal example

```r
library(sequencemap)

result <- run_sequence_map_workflow(
  fasta_path = "data/fasta/vp1_A_trimmed.fasta",
  metadata_path = "path/to/metadata.csv",
  boundary_path = "path/to/boundary.gpkg",
  output_dir = "outputs/example-run",
  sequence_id_col = "sequence_id",
  lon_col = "longitude",
  lat_col = "latitude",
  location_col = "location",
  crs_projected = "+proj=aea +lat_1=20 +lat_2=40 +lat_0=30 +lon_0=105 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
  grid_resolution = 10000,
  n_pca = 20,
  n_axes_model = 3,
  run_dapc = FALSE,
  overwrite = TRUE
)
```

## Tutorial and analysis workbooks

Render the general tutorial workbook from the repository root:

```r
quarto::quarto_render("analysis/sequence_map_workbook.qmd")
```

The tutorial uses bundled Vietnam example data and explains each major step. It
can render without INLA; mesh, model, and prediction sections are skipped when
INLA is unavailable.

The formal Vietnam analysis workbook is separate:

```r
quarto::quarto_render("analysis/vietnam_sequence_map_analysis.qmd")
```

That workbook is designed to run locally because fitting INLA models across all
FASTA alignments may take time.

## Output interpretation

PCA axis surfaces are predicted spatial patterns in genetic-position scores, not
direct estimates of nucleotide diversity. Location-level dispersion is an
ordination-space diversity summary and needs multiple sequences per location.
Predictions in unsampled areas are model-based interpolations rather than
observed genetic data.

GeoTIFFs and province summaries are reporting products derived from the
continuous prediction surface. `calculate_polygon_axis_metrics()` derives:

- `axis_centroid_distance`: predicted genetic-position displacement.
- `axis_surface_dispersion`: within-polygon predicted-axis heterogeneity.
- `mean_axis_uncertainty`: average posterior uncertainty.

These metrics should be interpreted alongside observed sequence support and
model uncertainty.

## Large output files

Generated results belong under `outputs/` and are ignored by Git. Avoid adding
large prediction grids, GeoPackages, model objects, rendered workbooks, or other
derived artifacts to the repository.

## Current limitations

- Gap and ambiguous sequence states are retained in QC summaries but are not
  modeled as separate nucleotide states by default.
- Spatial surfaces currently model one location-level ordination axis at a time.
- Ordination-space dispersion is not nucleotide diversity.
- PCA axes are alignment-specific unless a shared feature basis is implemented.
