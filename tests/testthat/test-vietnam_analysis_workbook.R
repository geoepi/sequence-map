test_that("the Vietnam analysis workbook is configured for safe local runs", {
  workbook_path <- testthat::test_path(
    "..", "..", "analysis", "vietnam_sequence_map_analysis.qmd"
  )

  expect_true(file.exists(workbook_path))
  workbook <- paste(readLines(workbook_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("run_models <- FALSE", workbook, fixed = TRUE))
  expect_true(grepl("here::here(", workbook, fixed = TRUE))
  expect_false(grepl("setwd(", workbook, fixed = TRUE))
  expect_true(grepl("list.files(", workbook, fixed = TRUE))
  expect_true(grepl("data/fasta", workbook, fixed = TRUE))
  expect_true(grepl("run_sequence_map_workflow(", workbook, fixed = TRUE))
  expect_true(grepl("write_axis_geotiffs(", workbook, fixed = TRUE))
  expect_true(grepl("aggregate_axis_rasters_to_polygons(", workbook, fixed = TRUE))
  expect_true(grepl("calculate_polygon_axis_metrics(", workbook, fixed = TRUE))
  expect_false(grepl("vp1_A_trimmed.fasta", workbook, fixed = TRUE))
  expect_true(grepl("not direct nucleotide diversity", workbook, fixed = TRUE))
})
