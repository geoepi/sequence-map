test_that("spatial preflight reports insufficient sampled locations", {
  expect_error(
    check_spatial_location_support(data.frame(location = c("a", "b"))),
    "at least 3 sampled locations"
  )
})

test_that("namespace exposes only the supported public API", {
  namespace_path <- testthat::test_path("..", "..", "NAMESPACE")
  namespace_lines <- readLines(namespace_path, warn = FALSE)
  exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", namespace_lines, value = TRUE))

  expected_exports <- c(
    "read_alignment",
    "read_metadata",
    "read_boundary",
    "validate_sequence_metadata",
    "validate_alignment_metadata_match",
    "alignment_to_variant_matrix",
    "run_sequence_pca",
    "run_sequence_dapc",
    "summarize_location_scores",
    "summarize_location_diversity",
    "prepare_spatial_points",
    "make_prediction_grid",
    "make_inla_mesh",
    "fit_axis_model",
    "fit_axis_models",
    "predict_axis_surface",
    "predict_axis_surfaces",
    "axis_surface_to_raster",
    "write_axis_geotiffs",
    "aggregate_axis_rasters_to_polygons",
    "calculate_polygon_axis_metrics",
    "run_sequence_map_workflow",
    "summarize_workflow_status",
    "write_workflow_status"
  )

  expect_setequal(exports, expected_exports)
  expect_false("hello_sequencemap" %in% exports)
})

test_that("package metadata and user documentation are release ready", {
  description_path <- testthat::test_path("..", "..", "DESCRIPTION")
  description <- readLines(description_path, warn = FALSE)
  readme_path <- testthat::test_path("..", "..", "README.md")
  readme <- paste(readLines(readme_path, warn = FALSE), collapse = "\n")

  expect_true(any(grepl("^Version: 0.1.0$", description)))
  expect_true(any(grepl("https://github.com/geoepi/sequence-map", description, fixed = TRUE)))
  expect_true(grepl("### Input contract", readme, fixed = TRUE))
  expect_true(grepl("workflow_status.csv", readme, fixed = TRUE))
  expect_true(grepl("images/sequence-map-sticker.png", readme, fixed = TRUE))
  expect_false(grepl("JMHumphreys/sequence-map", readme, fixed = TRUE))
})
