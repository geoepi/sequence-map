create_workflow_test_files <- function() {
  test_dir <- file.path(tempdir(), paste0("sequencemap-workflow-", as.integer(stats::runif(1, 1, 1e6))))
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)

  fasta_path <- file.path(test_dir, "alignment.fasta")
  writeLines(
    c(
      ">seq1", "ACGTAC",
      ">seq2", "ATGTAC",
      ">seq3", "GTATAC",
      ">seq4", "GTATTC"
    ),
    fasta_path
  )

  metadata_path <- file.path(test_dir, "metadata.csv")
  metadata <- data.frame(
    sequence_id = c("seq1", "seq2", "seq3", "seq4"),
    longitude = c(-93.2, -93.0, -92.8, -92.6),
    latitude = c(44.8, 44.9, 45.0, 45.1),
    stringsAsFactors = FALSE
  )
  utils::write.csv(metadata, metadata_path, row.names = FALSE)

  boundary <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(-93.5, 44.6),
        c(-92.3, 44.6),
        c(-92.3, 45.3),
        c(-93.5, 45.3),
        c(-93.5, 44.6)
      ))),
      crs = 4326
    )
  )
  boundary_path <- file.path(test_dir, "boundary.gpkg")
  sf::st_write(boundary, boundary_path, quiet = TRUE, delete_dsn = TRUE)

  list(
    test_dir = test_dir,
    fasta_path = fasta_path,
    metadata_path = metadata_path,
    boundary_path = boundary_path
  )
}

test_that("workflow runs PCA-first and returns expected outputs", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  files <- create_workflow_test_files()
  output_dir <- file.path(files$test_dir, "outputs")

  result <- run_sequence_map_workflow(
    fasta_path = files$fasta_path,
    metadata_path = files$metadata_path,
    boundary_path = files$boundary_path,
    output_dir = output_dir,
    crs_projected = 3857,
    grid_resolution = 20000,
    n_pca = 3,
    n_axes_model = 1,
    run_dapc = FALSE,
    verbose = FALSE
  )

  expect_s3_class(result, "sequencemap_workflow_result")
  expect_true(dir.exists(result$output_dir))
  expect_true(is.null(result$dapc))
  expect_true("PC1_mean_pred" %in% names(result$axis_surfaces$combined_surface))
  expect_true(is.list(result$workflow_status))

  required_files <- c(
    "alignment_site_summary.csv",
    "alignment_sequence_summary.csv",
    "pca_scores.csv",
    "pca_variance_summary.csv",
    "location_scores.csv",
    "location_diversity.csv",
    "axis_model_summaries.csv",
    "prediction_axis_surfaces.csv",
    "workflow_status.csv",
    "dapc_status.csv",
    "axis_model_status.csv",
    "axis_prediction_status.csv",
    "workflow_error_log.csv",
    "workflow_result.rds",
    "axis_models.rds",
    "axis_surfaces.rds"
  )
  expect_true(all(file.exists(file.path(result$output_dir, required_files))))
})

test_that("workflow does not require group_col when DAPC is disabled", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  files <- create_workflow_test_files()

  expect_s3_class(
    run_sequence_map_workflow(
      fasta_path = files$fasta_path,
      metadata_path = files$metadata_path,
      boundary_path = files$boundary_path,
      output_dir = file.path(files$test_dir, "outputs-no-dapc"),
      crs_projected = 3857,
      grid_resolution = 20000,
      n_pca = 3,
      n_axes_model = 1,
      run_dapc = FALSE,
      verbose = FALSE
    ),
    "sequencemap_workflow_result"
  )
})

test_that("workflow errors clearly when DAPC is requested without grouping information", {
  files <- create_workflow_test_files()

  expect_error(
    run_sequence_map_workflow(
      fasta_path = files$fasta_path,
      metadata_path = files$metadata_path,
      boundary_path = files$boundary_path,
      output_dir = file.path(files$test_dir, "outputs-dapc-error"),
      crs_projected = 3857,
      grid_resolution = 20000,
      run_dapc = TRUE,
      infer_dapc_groups = FALSE,
      verbose = FALSE
    ),
    "supply `group_col` or set `infer_dapc_groups = TRUE`"
  )
})

test_that("workflow creates a timestamped subdirectory when overwrite is false", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  files <- create_workflow_test_files()
  output_dir <- file.path(files$test_dir, "outputs-shared")

  first_run <- run_sequence_map_workflow(
    fasta_path = files$fasta_path,
    metadata_path = files$metadata_path,
    boundary_path = files$boundary_path,
    output_dir = output_dir,
    crs_projected = 3857,
    grid_resolution = 20000,
    n_pca = 3,
    n_axes_model = 1,
    verbose = FALSE
  )
  second_run <- run_sequence_map_workflow(
    fasta_path = files$fasta_path,
    metadata_path = files$metadata_path,
    boundary_path = files$boundary_path,
    output_dir = output_dir,
    crs_projected = 3857,
    grid_resolution = 20000,
    n_pca = 3,
    n_axes_model = 1,
    verbose = FALSE
  )

  expect_false(identical(first_run$output_dir, second_run$output_dir))
  expect_true(startsWith(second_run$output_dir, normalizePath(output_dir, winslash = "/", mustWork = TRUE)))
})
