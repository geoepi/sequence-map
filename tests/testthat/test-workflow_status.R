make_workflow_status_fixture <- function(include_errors = TRUE) {
  model_error <- if (include_errors) {
    list(
      response_col = "PC2_mean",
      message = "Model convergence failed.",
      class = c("simpleError", "error", "condition")
    )
  } else {
    NULL
  }
  prediction_error <- if (include_errors) {
    list(
      response_col = "PC2_mean",
      message = "Prediction matrix failed.",
      class = c("simpleError", "error", "condition")
    )
  } else {
    NULL
  }

  list(
    output_dir = file.path(tempdir(), "sequencemap-workflow-status"),
    alignment_qc = list(
      n_sequences = 8,
      n_sites = 120,
      n_retained_sites = 25
    ),
    location_scores = data.frame(location = c("a", "b")),
    prediction_grid = data.frame(grid_id = paste0("g", 1:4)),
    dapc = NULL,
    axis_models = list(
      response_cols = if (include_errors) c("PC1_mean", "PC2_mean") else "PC1_mean",
      models = list(
        PC1_mean = list(
          model_summary = list(
            n_observations = 8,
            n_mesh_vertices = 16,
            waic = 101,
            dic = 99
          )
        )
      ),
      errors = if (include_errors) list(PC2_mean = model_error) else list()
    ),
    axis_surfaces = list(
      surfaces = list(PC1_mean = data.frame(grid_id = paste0("g", 1:4))),
      errors = if (include_errors) list(PC2_mean = prediction_error) else list()
    ),
    settings = list(
      fasta_path = "alignment.fasta",
      run_dapc = FALSE
    )
  )
}

test_that("summarize_workflow_status records successful and failed axes", {
  status <- summarize_workflow_status(make_workflow_status_fixture())

  expect_equal(status$workflow_status$n_axis_models_requested, 2)
  expect_equal(status$workflow_status$n_axis_models_successful, 1)
  expect_equal(status$workflow_status$n_axis_models_failed, 1)
  expect_equal(status$workflow_status$n_axis_predictions_successful, 1)
  expect_equal(status$workflow_status$n_axis_predictions_failed, 1)
  expect_equal(status$dapc_status$status, "skipped")
  expect_equal(status$axis_model_status$status, c("success", "failed"))
  expect_equal(status$axis_prediction_status$status, c("success", "failed"))
  expect_true(all(c("axis_model", "axis_prediction") %in% status$error_log$stage))
  expect_true(all(c("stage", "response_col", "message", "class") %in% names(status$error_log)))
})

test_that("summarize_workflow_status marks skipped predictions", {
  workflow_result <- make_workflow_status_fixture(include_errors = FALSE)
  workflow_result$axis_models$response_cols <- c("PC1_mean", "PC2_mean")

  status <- summarize_workflow_status(workflow_result)

  expect_equal(status$axis_prediction_status$status, c("success", "skipped"))
  expect_equal(status$workflow_status$n_axis_predictions_failed, 0)
})

test_that("summarize_workflow_status records external stage errors", {
  status <- summarize_workflow_status(
    make_workflow_status_fixture(include_errors = FALSE),
    raster_outputs = simpleError("Raster export failed."),
    polygon_summaries = data.frame(polygon_id = "a")
  )

  expect_equal(status$workflow_status$raster_status, "failed")
  expect_equal(status$workflow_status$polygon_status, "success")
  expect_equal(status$error_log$stage, "raster")
})

test_that("write_workflow_status writes all tables including empty errors", {
  workflow_result <- make_workflow_status_fixture(include_errors = FALSE)
  status <- summarize_workflow_status(workflow_result)
  output_dir <- file.path(tempdir(), "sequencemap-workflow-status-files")
  files <- write_workflow_status(status, output_dir)

  expect_true(all(file.exists(files)))
  error_log <- utils::read.csv(files[["error_log"]], stringsAsFactors = FALSE)
  expect_equal(nrow(error_log), 0)
  expect_equal(
    names(error_log),
    c("stage", "response_col", "message", "class")
  )
})
