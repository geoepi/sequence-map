#' Run the end-to-end sequence-map workflow
#'
#' Runs the current PCA-first sequence-map pipeline from input files through
#' spatial axis-surface prediction. DAPC is optional and disabled by default.
#'
#' @param fasta_path Path to the FASTA alignment file.
#' @param metadata_path Path to the sequence metadata file.
#' @param boundary_path Path to a boundary readable by [sf::st_read()].
#' @param output_dir Output directory for workflow artifacts.
#' @param sequence_id_col Sequence identifier column name.
#' @param lon_col Longitude column name.
#' @param lat_col Latitude column name.
#' @param location_col Optional location identifier column.
#' @param group_col Optional grouping column for DAPC.
#' @param crs_input Input coordinate CRS.
#' @param crs_projected Projected CRS for spatial analysis.
#' @param grid_resolution Grid spacing in projected units.
#' @param n_pca Number of principal components to retain.
#' @param n_axes_model Number of PCA and optional DAPC axes to model.
#' @param run_dapc Logical; whether to run DAPC.
#' @param infer_dapc_groups Logical; whether to infer DAPC groups.
#' @param drop_invariant Logical; passed to [alignment_to_variant_matrix()].
#' @param iid_effects Optional iid-effect columns for axis models.
#' @param continue_on_error Logical; passed to batch model wrappers.
#' @param overwrite Logical; whether to reuse an existing output directory.
#' @param verbose Logical; whether to emit progress messages.
#'
#' @return A list of class `sequencemap_workflow_result`.
#' @export
run_sequence_map_workflow <- function(
  fasta_path,
  metadata_path,
  boundary_path,
  output_dir,
  sequence_id_col = "sequence_id",
  lon_col = "longitude",
  lat_col = "latitude",
  location_col = NULL,
  group_col = NULL,
  crs_input = 4326,
  crs_projected,
  grid_resolution,
  n_pca = 20,
  n_axes_model = 3,
  run_dapc = FALSE,
  infer_dapc_groups = FALSE,
  drop_invariant = TRUE,
  iid_effects = NULL,
  continue_on_error = TRUE,
  overwrite = FALSE,
  verbose = TRUE
) {
  if (missing(crs_projected) || is.null(crs_projected)) {
    stop("`crs_projected` is required for the workflow.", call. = FALSE)
  }
  if (missing(grid_resolution) || is.null(grid_resolution)) {
    stop("`grid_resolution` is required for the workflow.", call. = FALSE)
  }

  if (isTRUE(run_dapc) && is.null(group_col) && !isTRUE(infer_dapc_groups)) {
    stop(
      "When `run_dapc = TRUE`, supply `group_col` or set `infer_dapc_groups = TRUE`.",
      call. = FALSE
    )
  }

  realized_output_dir <- realize_workflow_output_dir(output_dir, overwrite = overwrite)
  workflow_message(verbose, sprintf("Writing workflow outputs to '%s'.", realized_output_dir))

  workflow_message(verbose, "Reading alignment, metadata, and boundary.")
  alignment <- read_alignment(fasta_path)
  metadata <- read_metadata(metadata_path)
  boundary <- read_boundary(boundary_path)

  workflow_message(verbose, "Validating sequence IDs and coordinates.")
  validate_alignment_metadata_match(alignment, metadata, id_col = sequence_id_col)
  validate_sequence_metadata(
    metadata = metadata,
    alignment = alignment,
    id_col = sequence_id_col,
    longitude_col = lon_col,
    latitude_col = lat_col
  )

  workflow_message(verbose, "Encoding alignment and running PCA.")
  alignment_dnabin <- normalize_workflow_alignment(alignment)
  variant_result <- alignment_to_variant_matrix(alignment_dnabin, drop_invariant = drop_invariant)
  pca_result <- run_sequence_pca(
    x = variant_result$variant_matrix,
    rank. = n_pca
  )

  dapc_result <- NULL
  if (isTRUE(run_dapc)) {
    workflow_message(verbose, "Running optional DAPC.")
    dapc_result <- run_sequence_dapc(
      x = variant_result$variant_matrix,
      metadata = metadata,
      id_col = sequence_id_col,
      grouping_col = group_col,
      infer_groups = infer_dapc_groups,
      n_pca = n_pca,
      n_da = n_axes_model
    )
  }

  workflow_message(verbose, "Combining sequence-level scores.")
  sequence_scores <- build_sequence_score_table(
    pca_result = pca_result,
    dapc_result = dapc_result
  )

  workflow_message(verbose, "Summarizing sampled locations.")
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

  modeled_response_cols <- select_modeled_response_cols(
    location_scores = location_scores,
    n_axes_model = n_axes_model,
    include_dapc = isTRUE(run_dapc)
  )

  workflow_message(verbose, "Preparing projected points, grid, and mesh.")
  points_sf <- prepare_spatial_points(
    data = location_scores,
    lon_col = lon_col,
    lat_col = lat_col,
    crs_input = crs_input,
    crs_projected = crs_projected,
    id_col = location_col,
    keep_cols = setdiff(names(location_scores), c(lon_col, lat_col))
  )
  prediction_grid <- make_prediction_grid(
    boundary = boundary,
    crs_projected = crs_projected,
    grid_resolution = grid_resolution,
    return = "sf"
  )
  check_spatial_workflow_dependencies()
  check_spatial_location_support(points_sf)
  mesh_result <- tryCatch(
    make_inla_mesh(
      points_sf = points_sf,
      boundary = boundary,
      prediction_grid = prediction_grid,
      max_edge = c(grid_resolution * 2, grid_resolution * 5),
      cutoff = grid_resolution / 5,
      offset = c(grid_resolution, grid_resolution * 2),
      crs_projected = crs_projected
    ),
    error = function(error) {
      stop(
        sprintf(
          "INLA mesh construction failed. Check the projected CRS, boundary, and sampled locations. Details: %s",
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )

  workflow_message(verbose, "Fitting batch axis models and predicting surfaces.")
  axis_models <- tryCatch(
    fit_axis_models(
      points_sf = points_sf,
      response_cols = modeled_response_cols,
      mesh = mesh_result,
      iid_effects = iid_effects,
      continue_on_error = continue_on_error,
      compute_criteria = TRUE,
      verbose = FALSE
    ),
    error = function(error) {
      stop(
        sprintf("INLA axis model fitting failed. Details: %s", conditionMessage(error)),
        call. = FALSE
      )
    }
  )
  axis_surfaces <- tryCatch(
    predict_axis_surfaces(
      axis_models = axis_models,
      prediction_grid = prediction_grid,
      continue_on_error = continue_on_error
    ),
    error = function(error) {
      stop(
        sprintf("INLA axis prediction failed. Details: %s", conditionMessage(error)),
        call. = FALSE
      )
    }
  )

  alignment_qc <- list(
    n_sequences = nrow(variant_result$sequence_summary),
    n_sites = nrow(variant_result$site_summary),
    n_retained_sites = sum(variant_result$site_summary$retained),
    n_variant_columns = ncol(variant_result$variant_matrix),
    drop_invariant = drop_invariant
  )

  workflow_result <- structure(
    list(
      output_dir = realized_output_dir,
      alignment_qc = alignment_qc,
      sequence_summary = variant_result$sequence_summary,
      site_summary = variant_result$site_summary,
      pca = pca_result,
      dapc = dapc_result,
      sequence_scores = sequence_scores,
      location_scores = location_scores,
      location_diversity = location_diversity,
      points_sf = points_sf,
      prediction_grid = prediction_grid,
      mesh = mesh_result,
      axis_models = axis_models,
      axis_surfaces = axis_surfaces,
      files_written = list(),
      settings = list(
        fasta_path = fasta_path,
        metadata_path = metadata_path,
        boundary_path = boundary_path,
        sequence_id_col = sequence_id_col,
        lon_col = lon_col,
        lat_col = lat_col,
        location_col = location_col,
        group_col = group_col,
        crs_input = crs_input,
        crs_projected = crs_projected,
        grid_resolution = grid_resolution,
        n_pca = n_pca,
        n_axes_model = n_axes_model,
        run_dapc = run_dapc,
        infer_dapc_groups = infer_dapc_groups,
        drop_invariant = drop_invariant,
        iid_effects = iid_effects,
        continue_on_error = continue_on_error,
        overwrite = overwrite
      )
    ),
    class = "sequencemap_workflow_result"
  )

  workflow_message(verbose, "Writing workflow outputs and diagnostics.")
  files_written <- write_workflow_outputs(
    workflow_result = workflow_result,
    output_dir = realized_output_dir,
    boundary = boundary,
    crs_input = crs_input,
    lon_col = lon_col,
    lat_col = lat_col,
    verbose = verbose
  )
  workflow_result$files_written <- files_written

  workflow_result$workflow_status <- summarize_workflow_status(workflow_result)
  status_files <- write_workflow_status(
    workflow_status = workflow_result$workflow_status,
    output_dir = realized_output_dir
  )
  workflow_result$files_written <- c(workflow_result$files_written, status_files)
  saveRDS(workflow_result, file.path(realized_output_dir, "workflow_result.rds"))

  workflow_result
}

#' Emit a workflow message when verbose is enabled
#'
#' @param verbose Logical flag.
#' @param text Message text.
#'
#' @return `NULL`, invisibly.
#' @keywords internal
workflow_message <- function(verbose, text) {
  if (isTRUE(verbose)) {
    message(text)
  }
  invisible(NULL)
}

#' Normalize workflow alignment input for encoding
#'
#' @param alignment Alignment object returned by [read_alignment()].
#'
#' @return A DNAbin alignment.
#' @keywords internal
normalize_workflow_alignment <- function(alignment) {
  if (inherits(alignment, "DNAbin")) {
    return(alignment)
  }

  normalized <- tryCatch(
    {
      dnabin <- ape::as.DNAbin(alignment)
      dnabin_matrix <- as.matrix(dnabin)
      class(dnabin_matrix) <- "DNAbin"
      dnabin_matrix
    },
    error = function(error) {
      stop(
        sprintf("Failed to coerce alignment to DNAbin for encoding: %s", error$message),
        call. = FALSE
      )
    }
  )

  normalized
}

#' Realize an output directory for a workflow run
#'
#' @param output_dir Requested output directory.
#' @param overwrite Logical flag.
#'
#' @return Character scalar path.
#' @keywords internal
realize_workflow_output_dir <- function(output_dir, overwrite = FALSE) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(output_dir, winslash = "/", mustWork = TRUE))
  }

  if (isTRUE(overwrite)) {
    return(normalizePath(output_dir, winslash = "/", mustWork = TRUE))
  }

  timestamp_dir <- file.path(output_dir, format(Sys.time(), "run_%Y%m%d_%H%M%S"))
  suffix <- 1L
  while (dir.exists(timestamp_dir)) {
    timestamp_dir <- file.path(output_dir, sprintf("%s_%02d", format(Sys.time(), "run_%Y%m%d_%H%M%S"), suffix))
    suffix <- suffix + 1L
  }
  dir.create(timestamp_dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(timestamp_dir, winslash = "/", mustWork = TRUE)
}

#' Build the sequence-level score table used downstream
#'
#' @param pca_result Result from [run_sequence_pca()].
#' @param dapc_result Optional result from [run_sequence_dapc()].
#'
#' @return Data frame.
#' @keywords internal
build_sequence_score_table <- function(pca_result, dapc_result = NULL) {
  pca_cols <- c("sequence_id", grep("^PC_?[0-9]+$", names(pca_result$scores), value = TRUE))
  sequence_scores <- pca_result$scores[, pca_cols, drop = FALSE]

  if (!is.null(dapc_result)) {
    dapc_cols <- c("sequence_id", grep("^LD_?[0-9]+$", names(dapc_result$scores), value = TRUE))
    dapc_scores <- dapc_result$scores[, dapc_cols, drop = FALSE]
    sequence_scores <- merge(sequence_scores, dapc_scores, by = "sequence_id", all = TRUE, sort = FALSE)
  }

  sequence_scores
}

#' Select modeled axis responses from location scores
#'
#' @param location_scores Location-level score summary table.
#' @param n_axes_model Number of axes per method.
#' @param include_dapc Logical; whether to include DAPC axes.
#'
#' @return Character vector.
#' @keywords internal
select_modeled_response_cols <- function(location_scores, n_axes_model = 3, include_dapc = FALSE) {
  pca_cols <- grep("^PC_?[0-9]+_mean$", names(location_scores), value = TRUE)
  pca_cols <- pca_cols[order(extract_axis_order(pca_cols))]
  selected <- utils::head(pca_cols, n_axes_model)

  if (isTRUE(include_dapc)) {
    dapc_cols <- grep("^LD_?[0-9]+_mean$", names(location_scores), value = TRUE)
    dapc_cols <- dapc_cols[order(extract_axis_order(dapc_cols))]
    selected <- c(selected, utils::head(dapc_cols, n_axes_model))
  }

  selected
}

#' Extract numeric axis order from column names
#'
#' @param x Character vector.
#'
#' @return Numeric vector.
#' @keywords internal
extract_axis_order <- function(x) {
  as.numeric(gsub("[^0-9]", "", x))
}

#' Write workflow CSV, RDS, and plot outputs
#'
#' @param workflow_result Workflow result object.
#' @param output_dir Output directory.
#' @param boundary Boundary sf object.
#' @param crs_input Input CRS.
#' @param lon_col Longitude column name.
#' @param lat_col Latitude column name.
#' @param verbose Logical.
#'
#' @return Named list of written files.
#' @keywords internal
write_workflow_outputs <- function(
  workflow_result,
  output_dir,
  boundary,
  crs_input,
  lon_col,
  lat_col,
  verbose = TRUE
) {
  files_written <- list()

  write_csv_output <- function(object, filename) {
    path <- file.path(output_dir, filename)
    utils::write.csv(object, path, row.names = FALSE)
    files_written[[filename]] <<- path
    workflow_message(verbose, sprintf("Wrote %s", filename))
  }

  save_rds_output <- function(object, filename) {
    path <- file.path(output_dir, filename)
    saveRDS(object, path)
    files_written[[filename]] <<- path
    workflow_message(verbose, sprintf("Wrote %s", filename))
  }

  write_csv_output(workflow_result$site_summary, "alignment_site_summary.csv")
  write_csv_output(workflow_result$sequence_summary, "alignment_sequence_summary.csv")
  write_csv_output(workflow_result$pca$scores, "pca_scores.csv")
  write_csv_output(workflow_result$pca$variance_summary, "pca_variance_summary.csv")

  if (!is.null(workflow_result$dapc)) {
    write_csv_output(workflow_result$dapc$scores, "dapc_scores.csv")
    write_csv_output(workflow_result$dapc$eigen_summary, "dapc_eigen_summary.csv")
  }

  write_csv_output(workflow_result$location_scores, "location_scores.csv")
  write_csv_output(workflow_result$location_diversity, "location_diversity.csv")
  write_csv_output(workflow_result$axis_models$model_summaries, "axis_model_summaries.csv")

  prediction_surface_csv <- prepare_prediction_surface_table(
    combined_surface = workflow_result$axis_surfaces$combined_surface,
    crs_input = crs_input,
    lon_col = lon_col,
    lat_col = lat_col
  )
  write_csv_output(prediction_surface_csv, "prediction_axis_surfaces.csv")

  if (!is.null(workflow_result$axis_surfaces$combined_surface)) {
    gpkg_path <- file.path(output_dir, "prediction_axis_surfaces.gpkg")
    sf::st_write(workflow_result$axis_surfaces$combined_surface, gpkg_path, delete_dsn = TRUE, quiet = TRUE)
    files_written[["prediction_axis_surfaces.gpkg"]] <- gpkg_path
  }

  save_rds_output(workflow_result$axis_models, "axis_models.rds")
  save_rds_output(workflow_result$axis_surfaces, "axis_surfaces.rds")

  workflow_result_copy <- workflow_result
  workflow_result_copy$files_written <- NULL
  save_rds_output(workflow_result_copy, "workflow_result.rds")

  maybe_write_pca_plot(workflow_result$pca$scores, output_dir)
  if (!is.null(workflow_result$dapc)) {
    maybe_write_dapc_plot(workflow_result$dapc$scores, output_dir)
  }
  maybe_write_mesh_plot(workflow_result$mesh, workflow_result$points_sf, output_dir)
  maybe_write_first_axis_map(
    surface = workflow_result$axis_surfaces$combined_surface,
    boundary = boundary,
    output_dir = output_dir
  )

  files_written
}

#' Prepare prediction surfaces for CSV output
#'
#' @param combined_surface Combined prediction sf object.
#' @param crs_input Input CRS.
#' @param lon_col Longitude column name.
#' @param lat_col Latitude column name.
#'
#' @return Data frame.
#' @keywords internal
prepare_prediction_surface_table <- function(combined_surface, crs_input, lon_col, lat_col) {
  if (is.null(combined_surface)) {
    return(data.frame())
  }

  prediction_table <- sf::st_drop_geometry(combined_surface)
  geographic_surface <- tryCatch(
    sf::st_transform(combined_surface, crs_input),
    error = function(error) NULL
  )

  if (!is.null(geographic_surface) && isTRUE(sf::st_is_longlat(geographic_surface))) {
    geographic_coords <- sf::st_coordinates(geographic_surface)
    prediction_table[[lon_col]] <- geographic_coords[, 1]
    prediction_table[[lat_col]] <- geographic_coords[, 2]
  }

  prediction_table
}

#' Write a simple PCA scatter plot when available
#'
#' @param pca_scores PCA score table.
#' @param output_dir Output directory.
#'
#' @return `NULL`, invisibly.
#' @keywords internal
maybe_write_pca_plot <- function(pca_scores, output_dir) {
  if (!all(c("PC1", "PC2") %in% names(pca_scores))) {
    return(invisible(NULL))
  }

  plot_data <- pca_scores
  plot_path <- file.path(output_dir, "pca_scatter.png")
  ggplot2::ggsave(
    filename = plot_path,
    plot = ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::geom_point() +
      ggplot2::theme_minimal(),
    width = 6,
    height = 5,
    dpi = 150
  )
  invisible(NULL)
}

#' Write a simple DAPC scatter plot when available
#'
#' @param dapc_scores DAPC score table.
#' @param output_dir Output directory.
#'
#' @return `NULL`, invisibly.
#' @keywords internal
maybe_write_dapc_plot <- function(dapc_scores, output_dir) {
  if (!all(c("LD1", "LD2") %in% names(dapc_scores))) {
    return(invisible(NULL))
  }

  plot_path <- file.path(output_dir, "dapc_scatter.png")
  ggplot2::ggsave(
    filename = plot_path,
    plot = ggplot2::ggplot(dapc_scores, ggplot2::aes(x = LD1, y = LD2)) +
      ggplot2::geom_point() +
      ggplot2::theme_minimal(),
    width = 6,
    height = 5,
    dpi = 150
  )
  invisible(NULL)
}

#' Write a mesh diagnostic plot if possible
#'
#' @param mesh_result Mesh list returned by [make_inla_mesh()].
#' @param points_sf Projected sf points.
#' @param output_dir Output directory.
#'
#' @return `NULL`, invisibly.
#' @keywords internal
maybe_write_mesh_plot <- function(mesh_result, points_sf, output_dir) {
  plot_path <- file.path(output_dir, "mesh_diagnostic.png")
  grDevices::png(plot_path, width = 1200, height = 900, res = 150)
  try({
    graphics::plot(mesh_result$mesh, asp = 1, main = "INLA mesh diagnostic")
    graphics::points(sf::st_coordinates(points_sf), col = "red", pch = 16)
  }, silent = TRUE)
  grDevices::dev.off()
  invisible(NULL)
}

#' Write a simple first-axis prediction map if available
#'
#' @param surface Combined prediction surface.
#' @param boundary Boundary sf object.
#' @param output_dir Output directory.
#'
#' @return `NULL`, invisibly.
#' @keywords internal
maybe_write_first_axis_map <- function(surface, boundary, output_dir) {
  if (is.null(surface)) {
    return(invisible(NULL))
  }

  prediction_cols <- grep("_pred$", names(surface), value = TRUE)
  if (length(prediction_cols) == 0L) {
    return(invisible(NULL))
  }

  plot_surface <- surface
  plot_surface$plot_value <- plot_surface[[prediction_cols[1]]]
  boundary_sf <- normalize_boundary(boundary)
  if (sf::st_crs(boundary_sf) != sf::st_crs(plot_surface)) {
    boundary_sf <- sf::st_transform(boundary_sf, sf::st_crs(plot_surface))
  }

  plot_path <- file.path(output_dir, "first_axis_map.png")
  ggplot2::ggsave(
    filename = plot_path,
    plot = ggplot2::ggplot() +
      ggplot2::geom_sf(data = plot_surface, ggplot2::aes(color = plot_value)) +
      ggplot2::geom_sf(data = boundary_sf, fill = NA, color = "black") +
      ggplot2::theme_minimal(),
    width = 6,
    height = 5,
    dpi = 150
  )
  invisible(NULL)
}
