test_that("tutorial script exists and references the expected workflow steps", {
  script_path <- testthat::test_path("..", "..", "scripts", "tutorial_sequence_map_example.R")

  expect_true(file.exists(script_path))

  script_lines <- readLines(script_path, warn = FALSE)
  script_text <- paste(script_lines, collapse = "\n")

  expect_true(any(grepl("^# ---- Step 1: Read inputs ----$", script_lines)))
  expect_true(any(grepl("^# ---- Step 10: INLA availability check ----$", script_lines)))
  expect_true(any(grepl("^# ---- Step 14: Compare to full wrapper ----$", script_lines)))

  expected_functions <- c(
    "read_alignment\\(",
    "read_metadata\\(",
    "read_boundary\\(",
    "validate_sequence_metadata\\(",
    "validate_alignment_metadata_match\\(",
    "alignment_to_variant_matrix\\(",
    "run_sequence_pca\\(",
    "run_sequence_dapc\\(",
    "summarize_location_scores\\(",
    "summarize_location_diversity\\(",
    "prepare_spatial_points\\(",
    "make_prediction_grid\\(",
    "make_inla_mesh\\(",
    "fit_axis_models\\(",
    "predict_axis_surfaces\\(",
    "run_sequence_map_workflow\\("
  )

  for (pattern in expected_functions) {
    expect_true(grepl(pattern, script_text))
  }

  expect_false(any(grepl("^source\\(", script_lines)))
})
