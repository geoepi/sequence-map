test_that("the tutorial workbook contains the sequence-map workflow", {
  workbook_path <- testthat::test_path("..", "..", "analysis", "sequence_map_workbook.qmd")

  expect_true(file.exists(workbook_path))

  workbook <- paste(readLines(workbook_path, warn = FALSE), collapse = "\n")
  expected_sections <- c(
    "## Interpretation first",
    "## Setup",
    "## Read and validate inputs",
    "## Alignment encoding diagnostics",
    "## Principal component analysis",
    "## Optional discriminant analysis of principal components",
    "## Location summaries",
    "## Spatial preparation",
    "## INLA availability",
    "## Mesh diagnostics",
    "## Axis models",
    "## Spatial predictions",
    "## Reproducibility"
  )

  for (section in expected_sections) {
    expect_true(grepl(section, workbook, fixed = TRUE))
  }

  expected_functions <- c(
    "read_alignment(",
    "alignment_to_variant_matrix(",
    "run_sequence_pca(",
    "summarize_location_scores(",
    "prepare_spatial_points(",
    "make_prediction_grid(",
    "make_inla_mesh(",
    "fit_axis_models(",
    "predict_axis_surfaces("
  )

  for (function_name in expected_functions) {
    expect_true(grepl(function_name, workbook, fixed = TRUE))
  }

  expect_true(grepl("write_outputs <- FALSE", workbook, fixed = TRUE))
  expect_true(grepl("requireNamespace(\"INLA\", quietly = TRUE)", workbook, fixed = TRUE))
  expect_true(grepl("eval: !expr inla_available", workbook, fixed = TRUE))
})
