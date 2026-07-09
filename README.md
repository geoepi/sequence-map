<!--
Add the package sticker at images/sequence-map-sticker.png, then uncomment:
<p align="center">
  <img src="images/sequence-map-sticker.png" width="220" alt="sequence-map hex sticker">
</p>
-->

# sequencemap

## What sequence-map does

`sequencemap` converts aligned viral sequences and sampling coordinates into
PCA-first, location-level summaries and optional spatial predictions of
ordination-axis scores. It supports alignment quality control, PCA, optional
DAPC, location aggregation, projected prediction grids, and Gaussian
SPDE/INLA axis models.

## Inputs

- A FASTA alignment with unique sequence identifiers.
- Sequence metadata with matching identifiers and numeric longitude/latitude.
- A boundary readable by `sf` for prediction-grid construction.

## Core workflow

`read_alignment()` and input validation -> `alignment_to_variant_matrix()` ->
`run_sequence_pca()` -> `summarize_location_scores()` ->
`prepare_spatial_points()` -> `make_prediction_grid()` -> optional
`make_inla_mesh()`, `fit_axis_models()`, and `predict_axis_surfaces()` ->
`axis_surface_to_raster()` -> `aggregate_axis_rasters_to_polygons()`.

DAPC is opt-in because it requires a biologically justified grouping decision
or explicitly exploratory inferred groups.

## Tutorial workbook

Render `analysis/sequence_map_workbook.qmd` from the repository root for a
diagnostic walkthrough using the bundled Vietnam example data. The workbook
does not write files by default and can render without INLA; mesh, model, and
prediction sections are skipped when INLA is unavailable.

## Output interpretation

PCA axis surfaces are predicted spatial patterns in genetic-position scores,
not direct estimates of nucleotide diversity. Location-level dispersion is an
ordination-space diversity summary and needs multiple sequences per location.
Predictions in unsampled areas are model-based interpolations rather than
observed genetic data. GeoTIFFs and province summaries are reporting products
of the continuous surface, not separate province-level models or nucleotide
diversity estimates. `calculate_polygon_axis_metrics()` derives
`axis_centroid_distance` (predicted genetic-position displacement),
`axis_surface_dispersion` (within-polygon predicted-axis heterogeneity), and
`mean_axis_uncertainty` (average posterior uncertainty) from those summaries.

## INLA note

INLA is listed under `Suggests` because its installation is nonstandard. Only
the mesh and spatial-model functions require it. All spatial modeling uses a
projected CRS in metric units.

## Large output files

Generated results belong under `outputs/` and are ignored by Git. Avoid adding
large prediction grids, GeoPackages, model objects, rendered workbooks, or
other derived artifacts to the repository.

## Current limitations

- The default encoder represents gap and ambiguous states as zeros while
  retaining their QC counts.
- Spatial surfaces currently model one location-level ordination axis at a
  time; they do not fit raw sequence-level scores by default.
- Ordination-space dispersion is not nucleotide diversity.
- The current package does not yet provide multiaxis or diversity-surface
  models.
