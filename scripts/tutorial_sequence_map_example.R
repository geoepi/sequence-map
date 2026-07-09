# Tutorial: sequence-map step-by-step example
# Expected inputs:
# - a FASTA alignment
# - sequence metadata with IDs and coordinates
# - a spatial boundary readable by sf
# Expected outputs:
# - intermediate CSV summaries, diagnostic plots, and optional spatial predictions
# Note:
# - INLA is only required for the mesh and spatial modeling sections
# - this script is designed to be run section-by-section interactively

# ---- User settings ----

# If repository example files exist, use them by default.
repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
example_fasta_path <- file.path(repo_root, "data", "fasta", "vp1_A_trimmed.fasta")
example_metadata_rda <- file.path(repo_root, "data", "metadata", "prep-for-seqs.Rda")
example_boundary_rda <- file.path(repo_root, "data", "spatial", "vn-spat.Rda")

fasta_path <- if (file.exists(example_fasta_path)) example_fasta_path else "path/to/alignment.fasta"
metadata_path <- if (file.exists(example_metadata_rda)) example_metadata_rda else "path/to/metadata.csv"
boundary_path <- if (file.exists(example_boundary_rda)) example_boundary_rda else "path/to/boundary.gpkg"
output_dir <- file.path(repo_root, "outputs", "tutorial-sequence-map")

sequence_id_col <- "sequence_id"
lon_col <- "longitude"
lat_col <- "latitude"
location_col <- "location"
group_col <- NULL
crs_input <- 4326
crs_projected <- 8859
grid_resolution <- 10000
n_pca <- 20
n_axes_model <- 3
run_dapc <- FALSE
infer_dapc_groups <- FALSE
drop_invariant <- TRUE

# ---- Package loading and output directory setup ----

library(sequencemap)
library(dplyr)
library(ggplot2)
library(sf)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Repository metadata and boundary examples are currently stored as .Rda files.
# Convert them to temporary tutorial files so the exported readers can be used.
materialize_rda_example <- function(path, object_class = c("data.frame", "sf")) {
  object_class <- match.arg(object_class)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  objects <- mget(ls(env), envir = env)

  if (object_class == "data.frame") {
    candidate_name <- names(objects)[vapply(objects, function(x) inherits(x, "data.frame"), logical(1))][1]
    if (is.na(candidate_name)) {
      stop("No data.frame object was found in metadata .Rda example.")
    }
    tmp_path <- tempfile(fileext = ".csv")
    utils::write.csv(objects[[candidate_name]], tmp_path, row.names = FALSE)
    return(tmp_path)
  }

  candidate_name <- names(objects)[vapply(objects, function(x) inherits(x, "sf"), logical(1))][1]
  if (is.na(candidate_name)) {
    stop("No sf object was found in boundary .Rda example.")
  }
  tmp_path <- tempfile(fileext = ".gpkg")
  sf::st_write(objects[[candidate_name]], tmp_path, quiet = TRUE, delete_dsn = TRUE)
  tmp_path
}

materialize_repository_metadata <- function(metadata_rda, boundary_rda) {
  metadata_env <- new.env(parent = emptyenv())
  boundary_env <- new.env(parent = emptyenv())
  load(metadata_rda, envir = metadata_env)
  load(boundary_rda, envir = boundary_env)

  if (!exists("seq_meta_reg", envir = metadata_env, inherits = FALSE)) {
    stop("Repository metadata example must contain `seq_meta_reg`.")
  }
  if (!exists("vn_provinces", envir = boundary_env, inherits = FALSE)) {
    stop("Repository boundary example must contain `vn_provinces`.")
  }

  metadata <- get("seq_meta_reg", envir = metadata_env)
  provinces <- get("vn_provinces", envir = boundary_env)
  province_points <- sf::st_transform(
    sf::st_point_on_surface(sf::st_geometry(sf::st_transform(provinces, 3857))),
    4326
  )
  coordinates <- sf::st_coordinates(province_points)
  location_coordinates <- data.frame(
    location = provinces$prov_eng,
    longitude = coordinates[, "X"],
    latitude = coordinates[, "Y"],
    stringsAsFactors = FALSE
  )

  match_index <- match(metadata$location, location_coordinates$location)
  metadata$longitude <- location_coordinates$longitude[match_index]
  metadata$latitude <- location_coordinates$latitude[match_index]
  names(metadata)[names(metadata) == "accession"] <- "sequence_id"

  tmp_path <- tempfile(fileext = ".csv")
  utils::write.csv(metadata, tmp_path, row.names = FALSE)
  tmp_path
}

if (tolower(tools::file_ext(metadata_path)) == "rda") {
  metadata_path <- materialize_repository_metadata(metadata_path, boundary_path)
}
if (tolower(tools::file_ext(boundary_path)) == "rda") {
  boundary_path <- materialize_rda_example(boundary_path, object_class = "sf")
}

# ---- Step 1: Read inputs ----

# Read the alignment, metadata, and boundary into memory.
alignment <- read_alignment(fasta_path)
metadata <- read_metadata(metadata_path)
boundary <- read_boundary(boundary_path)

alignment_dnabin <- alignment
character_alignment <- ape::as.character.DNAbin(alignment_dnabin)
if (is.list(character_alignment)) {
  character_alignment <- do.call(rbind, character_alignment)
}
alignment_dimensions <- dim(character_alignment)
cat("Alignment dimensions:", alignment_dimensions[1], "sequences x", alignment_dimensions[2], "sites\n")
cat("Metadata dimensions:", nrow(metadata), "rows x", ncol(metadata), "columns\n")
cat("Boundary CRS:", sf::st_crs(boundary)$input, "\n")
print(sf::st_bbox(boundary))

# ---- Step 2: Validate inputs ----

# Confirm sequence IDs align and coordinate fields are usable.
validate_sequence_metadata(
  metadata = metadata,
  alignment = alignment,
  id_col = sequence_id_col,
  longitude_col = lon_col,
  latitude_col = lat_col
)
validate_alignment_metadata_match(alignment, metadata, id_col = sequence_id_col)
message("Input validation completed successfully.")

# ---- Step 3: Convert alignment to variant matrix ----

# Encode canonical sequence variation for ordination.
variant_result <- alignment_to_variant_matrix(
  alignment = alignment_dnabin,
  drop_invariant = drop_invariant
)
variant_matrix <- variant_result$variant_matrix
site_summary <- variant_result$site_summary
sequence_summary <- variant_result$sequence_summary

cat("Number of sequences:", nrow(variant_matrix), "\n")
cat("Number of retained sites:", sum(site_summary$retained), "\n")
cat("Number of encoded columns:", ncol(variant_matrix), "\n")
print(table(site_summary$site_class))

utils::write.csv(site_summary, file.path(output_dir, "alignment_site_summary.csv"), row.names = FALSE)
utils::write.csv(sequence_summary, file.path(output_dir, "alignment_sequence_summary.csv"), row.names = FALSE)

# ---- Step 4: Run PCA ----

# Run PCA on the encoded alignment matrix.
pca_result <- run_sequence_pca(variant_matrix, rank. = n_pca)
print(utils::head(pca_result$scores))
print(utils::head(pca_result$variance_summary, 10))

utils::write.csv(pca_result$scores, file.path(output_dir, "pca_scores.csv"), row.names = FALSE)
utils::write.csv(pca_result$variance_summary, file.path(output_dir, "pca_variance_summary.csv"), row.names = FALSE)

if (all(c("PC1", "PC2") %in% names(pca_result$scores))) {
  ggplot2::ggsave(
    filename = file.path(output_dir, "pca_scatter.png"),
    plot = ggplot2::ggplot(pca_result$scores, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::geom_point() +
      ggplot2::theme_minimal(),
    width = 6,
    height = 5,
    dpi = 150
  )
}

# ---- Step 5: Optional DAPC ----

dapc_result <- NULL
if (isTRUE(run_dapc)) {
  if (is.null(group_col) && !isTRUE(infer_dapc_groups)) {
    stop("Set `group_col` or `infer_dapc_groups = TRUE` before running DAPC.")
  }

  dapc_result <- run_sequence_dapc(
    x = variant_matrix,
    metadata = metadata,
    id_col = sequence_id_col,
    grouping_col = group_col,
    infer_groups = infer_dapc_groups,
    n_pca = n_pca,
    n_da = n_axes_model
  )

  utils::write.csv(dapc_result$scores, file.path(output_dir, "dapc_scores.csv"), row.names = FALSE)
  utils::write.csv(dapc_result$eigen_summary, file.path(output_dir, "dapc_eigen_summary.csv"), row.names = FALSE)

  if (all(c("LD1", "LD2") %in% names(dapc_result$scores))) {
    ggplot2::ggsave(
      filename = file.path(output_dir, "dapc_scatter.png"),
      plot = ggplot2::ggplot(dapc_result$scores, ggplot2::aes(x = LD1, y = LD2)) +
        ggplot2::geom_point() +
        ggplot2::theme_minimal(),
      width = 6,
      height = 5,
      dpi = 150
    )
  }
} else {
  message("DAPC was skipped. Set `run_dapc <- TRUE` to include it.")
}

# ---- Step 6: Combine sequence score tables ----

# Merge PCA scores with optional DAPC scores, keeping one row per sequence.
sequence_scores <- pca_result$scores %>%
  dplyr::select(sequence_id, dplyr::matches("^PC_?[0-9]+$"))

if (!is.null(dapc_result)) {
  sequence_scores <- sequence_scores %>%
    dplyr::left_join(
      dapc_result$scores %>% dplyr::select(sequence_id, dplyr::matches("^LD_?[0-9]+$")),
      by = "sequence_id"
    )
}

print(names(sequence_scores))
cat("Sequences represented in score table:", nrow(sequence_scores), "\n")

# ---- Step 7: Summarize sampled locations ----

# Aggregate sequence-level scores to sampled locations.
location_scores <- summarize_location_scores(
  scores = sequence_scores,
  metadata = metadata,
  sequence_id_col = sequence_id_col,
  location_col = location_col,
  lon_col = lon_col,
  lat_col = lat_col
)
location_diversity <- summarize_location_diversity(
  scores = sequence_scores,
  metadata = metadata,
  sequence_id_col = sequence_id_col,
  location_col = location_col,
  lon_col = lon_col,
  lat_col = lat_col
)

cat("Number of sampled locations:", nrow(location_scores), "\n")
cat("Locations with >= 2 sequences:", sum(location_diversity$n_sequences >= 2, na.rm = TRUE), "\n")
print(utils::head(location_scores))
print(utils::head(location_diversity))

utils::write.csv(location_scores, file.path(output_dir, "location_scores.csv"), row.names = FALSE)
utils::write.csv(location_diversity, file.path(output_dir, "location_diversity.csv"), row.names = FALSE)

# ---- Step 8: Prepare projected spatial points ----

# Convert sampled locations into projected sf points for mesh/model construction.
points_sf <- prepare_spatial_points(
  data = location_scores,
  lon_col = lon_col,
  lat_col = lat_col,
  crs_input = crs_input,
  crs_projected = crs_projected,
  keep_cols = setdiff(names(location_scores), c(lon_col, lat_col))
)

cat("Projected CRS:", sf::st_crs(points_sf)$input, "\n")
cat("X range:", paste(range(points_sf$x), collapse = " to "), "\n")
cat("Y range:", paste(range(points_sf$y), collapse = " to "), "\n")
cat("Number of sampled points:", nrow(points_sf), "\n")

# ---- Step 9: Create prediction grid ----

# Build a regular grid of prediction points inside the boundary.
prediction_grid <- make_prediction_grid(
  boundary = boundary,
  crs_projected = crs_projected,
  grid_resolution = grid_resolution,
  return = "sf"
)

cat("Number of grid points:", nrow(prediction_grid), "\n")
cat("Grid resolution:", grid_resolution, "\n")
cat("Grid CRS:", sf::st_crs(prediction_grid)$input, "\n")

sf::st_write(prediction_grid, file.path(output_dir, "prediction_grid.gpkg"), quiet = TRUE, delete_dsn = TRUE)

# ---- Step 10: INLA availability check ----

if (!requireNamespace("INLA", quietly = TRUE)) {
  message("INLA is not installed. Stopping before mesh/model fitting.")
  message("Install INLA locally, then rerun from Step 10.")
  stop("INLA unavailable; tutorial completed through non-INLA steps.")
}

# ---- Step 11: Create INLA mesh ----

# Build the SPDE mesh over sampled points and prediction domain.
mesh_result <- make_inla_mesh(
  points_sf = points_sf,
  boundary = boundary,
  prediction_grid = prediction_grid,
  max_edge = c(grid_resolution * 2, grid_resolution * 5),
  cutoff = grid_resolution / 5,
  offset = c(grid_resolution, grid_resolution * 2),
  crs_projected = crs_projected
)

print(mesh_result$mesh_summary)
cat("Number of mesh vertices:", mesh_result$mesh_summary$n_vertices, "\n")
cat("Number of mesh triangles:", mesh_result$mesh_summary$n_triangles, "\n")

grDevices::png(file.path(output_dir, "mesh_diagnostic.png"), width = 1200, height = 900, res = 150)
plot(mesh_result$mesh, asp = 1, main = "INLA mesh diagnostic")
points(sf::st_coordinates(points_sf), col = "red", pch = 16)
grDevices::dev.off()

# ---- Step 12: Fit batch axis models ----

# Select the first requested PCA axes and optional DAPC axes for spatial modeling.
response_cols <- names(location_scores)[grepl("^(PC|PC_|LD|LD_).+_mean$", names(location_scores))]
response_cols <- response_cols[order(as.numeric(gsub("[^0-9]", "", response_cols)))]
response_cols <- utils::head(response_cols, n_axes_model)

axis_models <- fit_axis_models(
  points_sf = points_sf,
  response_cols = response_cols,
  mesh = mesh_result,
  continue_on_error = TRUE
)

print(axis_models$successful_response_cols)
print(axis_models$failed_response_cols)
print(axis_models$model_summaries)

utils::write.csv(axis_models$model_summaries, file.path(output_dir, "axis_model_summaries.csv"), row.names = FALSE)
saveRDS(axis_models, file.path(output_dir, "axis_models.rds"))

# ---- Step 13: Predict axis surfaces ----

# Project each fitted axis model to the common prediction grid.
axis_surfaces <- predict_axis_surfaces(
  axis_models = axis_models,
  prediction_grid = prediction_grid,
  continue_on_error = TRUE
)

print(names(axis_surfaces$surfaces))
print(names(axis_surfaces$combined_surface))

prediction_surface_table <- sf::st_drop_geometry(axis_surfaces$combined_surface)
utils::write.csv(prediction_surface_table, file.path(output_dir, "prediction_axis_surfaces.csv"), row.names = FALSE)
sf::st_write(axis_surfaces$combined_surface, file.path(output_dir, "prediction_axis_surfaces.gpkg"), quiet = TRUE, delete_dsn = TRUE)
saveRDS(axis_surfaces, file.path(output_dir, "axis_surfaces.rds"))

first_prediction_col <- grep("_pred$", names(axis_surfaces$combined_surface), value = TRUE)[1]
if (!is.na(first_prediction_col)) {
  plot_surface <- axis_surfaces$combined_surface
  plot_surface$plot_value <- plot_surface[[first_prediction_col]]
  boundary_projected <- sf::st_transform(boundary, sf::st_crs(plot_surface))

  ggplot2::ggsave(
    filename = file.path(output_dir, "first_axis_map.png"),
    plot = ggplot2::ggplot() +
      ggplot2::geom_sf(data = plot_surface, ggplot2::aes(color = plot_value)) +
      ggplot2::geom_sf(data = boundary_projected, fill = NA, color = "black") +
      ggplot2::theme_minimal(),
    width = 6,
    height = 5,
    dpi = 150
  )
}

# ---- Step 14: Compare to full wrapper ----

# The full pipeline can also be run with a single wrapper call:
#
# workflow_result <- run_sequence_map_workflow(
#   fasta_path = fasta_path,
#   metadata_path = metadata_path,
#   boundary_path = boundary_path,
#   output_dir = output_dir,
#   sequence_id_col = sequence_id_col,
#   lon_col = lon_col,
#   lat_col = lat_col,
#   location_col = location_col,
#   group_col = group_col,
#   crs_input = crs_input,
#   crs_projected = crs_projected,
#   grid_resolution = grid_resolution,
#   n_pca = n_pca,
#   n_axes_model = n_axes_model,
#   run_dapc = run_dapc,
#   infer_dapc_groups = infer_dapc_groups,
#   drop_invariant = drop_invariant
# )
