#' Summarize workflow execution status
#'
#' Builds human-readable audit tables from a sequence-map workflow result.
#'
#' @param workflow_result A sequencemap workflow result.
#' @param raster_outputs Optional raster export result or error.
#' @param polygon_summaries Optional polygon summary result or error.
#' @param polygon_metrics Optional polygon metric result or error.
#'
#' @return A named list of audit data frames.
#' @export
summarize_workflow_status <- function(
  workflow_result,
  raster_outputs = NULL,
  polygon_summaries = NULL,
  polygon_metrics = NULL
) {
  if (!is.list(workflow_result)) {
    stop("'workflow_result' must be a workflow result list.", call. = FALSE)
  }
  settings <- workflow_result$settings
  if (is.null(settings)) settings <- list()
  alignment_qc <- workflow_result$alignment_qc
  if (is.null(alignment_qc)) alignment_qc <- list()

  dapc_status <- build_dapc_status(workflow_result, settings)
  axis_model_status <- build_axis_model_status(workflow_result$axis_models)
  axis_prediction_status <- build_axis_prediction_status(
    workflow_result$axis_models,
    workflow_result$axis_surfaces
  )
  error_log <- combine_workflow_error_log(
    dapc_status,
    axis_model_status,
    axis_prediction_status,
    raster_outputs,
    polygon_summaries,
    polygon_metrics
  )

  workflow_status <- data.frame(
    output_dir = scalar_character(workflow_result$output_dir),
    fasta_path = scalar_character(settings$fasta_path),
    n_sequences = scalar_numeric(alignment_qc$n_sequences),
    n_sites = scalar_numeric(alignment_qc$n_sites),
    n_retained_sites = scalar_numeric(alignment_qc$n_retained_sites),
    n_sampled_locations = nrow_or_na(workflow_result$location_scores),
    n_prediction_grid_cells = nrow_or_na(workflow_result$prediction_grid),
    run_dapc = isTRUE(settings$run_dapc),
    dapc_status = dapc_status$status[1],
    n_axis_models_requested = nrow(axis_model_status),
    n_axis_models_successful = sum(axis_model_status$status == "success"),
    n_axis_models_failed = sum(axis_model_status$status == "failed"),
    n_axis_predictions_successful = sum(axis_prediction_status$status == "success"),
    n_axis_predictions_failed = sum(axis_prediction_status$status == "failed"),
    raster_status = external_stage_status(raster_outputs),
    polygon_status = external_stage_status(polygon_summaries),
    polygon_metrics_status = external_stage_status(polygon_metrics),
    stringsAsFactors = FALSE
  )

  list(
    workflow_status = workflow_status,
    dapc_status = dapc_status,
    axis_model_status = axis_model_status,
    axis_prediction_status = axis_prediction_status,
    error_log = error_log
  )
}

#' Write workflow execution status tables
#'
#' @param workflow_status Audit list returned by summarize_workflow_status.
#' @param output_dir Optional output directory override.
#'
#' @return A named list of CSV file paths.
#' @export
write_workflow_status <- function(workflow_status, output_dir = NULL) {
  table_names <- c(
    "workflow_status",
    "dapc_status",
    "axis_model_status",
    "axis_prediction_status",
    "error_log"
  )
  if (!is.list(workflow_status) || !all(table_names %in% names(workflow_status))) {
    stop("'workflow_status' must be returned by summarize_workflow_status.", call. = FALSE)
  }
  if (is.null(output_dir)) output_dir <- workflow_status$workflow_status$output_dir[1]
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(output_dir)) {
    stop("'output_dir' must be supplied or recorded in workflow_status.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  filenames <- c(
    workflow_status = "workflow_status.csv",
    dapc_status = "dapc_status.csv",
    axis_model_status = "axis_model_status.csv",
    axis_prediction_status = "axis_prediction_status.csv",
    error_log = "workflow_error_log.csv"
  )
  paths <- file.path(output_dir, unname(filenames))
  names(paths) <- names(filenames)
  for (table_name in names(filenames)) {
    utils::write.csv(workflow_status[[table_name]], paths[[table_name]], row.names = FALSE)
  }
  paths
}

#' Build DAPC audit status
#' @keywords internal
build_dapc_status <- function(workflow_result, settings) {
  result <- workflow_result$dapc
  requested <- isTRUE(settings$run_dapc)
  status <- if (!requested) "skipped" else if (inherits(result, "error") || is.null(result)) "failed" else "success"
  message <- if (inherits(result, "error")) {
    conditionMessage(result)
  } else if (requested && is.null(result)) {
    "DAPC was requested but no result was recorded."
  } else {
    NA_character_
  }
  data.frame(
    run_dapc = requested,
    status = status,
    n_sequences = if (!is.null(result$scores)) nrow(result$scores) else NA_integer_,
    n_axes = if (!is.null(result$eigen_summary)) nrow(result$eigen_summary) else NA_integer_,
    error_message = message,
    error_class = if (inherits(result, "error")) paste(class(result), collapse = "/") else NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Build per-axis model audit status
#' @keywords internal
build_axis_model_status <- function(axis_models) {
  requested <- if (!is.null(axis_models$response_cols)) axis_models$response_cols else character()
  if (length(requested) == 0L) return(empty_axis_model_status())
  do.call(rbind, lapply(requested, function(response_col) {
    model <- axis_models$models[[response_col]]
    error <- axis_models$errors[[response_col]]
    summary <- if (!is.null(model$model_summary)) model$model_summary else list()
    data.frame(
      response_col = response_col,
      status = if (!is.null(model)) "success" else if (!is.null(error)) "failed" else "skipped",
      n_observations = scalar_numeric(summary$n_observations),
      n_mesh_vertices = scalar_numeric(summary$n_mesh_vertices),
      waic = scalar_numeric(summary$waic),
      dic = scalar_numeric(summary$dic),
      error_message = if (!is.null(error$message)) as.character(error$message) else NA_character_,
      error_class = if (!is.null(error$class)) paste(error$class, collapse = "/") else NA_character_,
      stringsAsFactors = FALSE
    )
  }))
}

#' Build per-axis prediction audit status
#' @keywords internal
build_axis_prediction_status <- function(axis_models, axis_surfaces) {
  requested <- if (!is.null(axis_models$response_cols)) axis_models$response_cols else character()
  if (length(requested) == 0L) return(empty_axis_prediction_status())
  do.call(rbind, lapply(requested, function(response_col) {
    surface <- axis_surfaces$surfaces[[response_col]]
    error <- axis_surfaces$errors[[response_col]]
    model_error <- axis_models$errors[[response_col]]
    status <- if (!is.null(surface)) "success" else if (!is.null(error)) "failed" else "skipped"
    message <- if (!is.null(error$message)) {
      as.character(error$message)
    } else if (!is.null(model_error$message)) {
      paste0("Model fitting did not succeed: ", model_error$message)
    } else {
      NA_character_
    }
    data.frame(
      response_col = response_col,
      status = status,
      n_prediction_cells = nrow_or_na(surface),
      error_message = message,
      error_class = if (!is.null(error$class)) {
        paste(error$class, collapse = "/")
      } else if (!is.null(model_error$class)) {
        paste(model_error$class, collapse = "/")
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  }))
}

#' Combine recorded workflow errors
#' @keywords internal
combine_workflow_error_log <- function(
  dapc_status,
  axis_model_status,
  axis_prediction_status,
  raster_outputs,
  polygon_summaries,
  polygon_metrics
) {
  rows <- list()
  if (dapc_status$status[1] == "failed") {
    rows[[length(rows) + 1L]] <- workflow_error_row(
      "dapc", NA_character_, dapc_status$error_message[1], dapc_status$error_class[1]
    )
  }
  for (index in which(axis_model_status$status == "failed")) {
    rows[[length(rows) + 1L]] <- workflow_error_row(
      "axis_model",
      axis_model_status$response_col[index],
      axis_model_status$error_message[index],
      axis_model_status$error_class[index]
    )
  }
  for (index in which(axis_prediction_status$status == "failed")) {
    rows[[length(rows) + 1L]] <- workflow_error_row(
      "axis_prediction",
      axis_prediction_status$response_col[index],
      axis_prediction_status$error_message[index],
      axis_prediction_status$error_class[index]
    )
  }
  external_values <- list(
    list(stage = "raster", value = raster_outputs),
    list(stage = "polygon", value = polygon_summaries),
    list(stage = "polygon", value = polygon_metrics)
  )
  for (item in external_values) {
    value <- item$value
    if (inherits(value, "error")) {
      rows[[length(rows) + 1L]] <- workflow_error_row(
        item$stage,
        NA_character_,
        conditionMessage(value),
        paste(class(value), collapse = "/")
      )
    }
  }
  if (length(rows) == 0L) return(empty_workflow_error_log())
  do.call(rbind, rows)
}

#' Create one workflow error row
#' @keywords internal
workflow_error_row <- function(stage, response_col, message, class) {
  data.frame(
    stage = stage,
    response_col = response_col,
    message = message,
    class = class,
    stringsAsFactors = FALSE
  )
}

#' Get an external-stage status
#' @keywords internal
external_stage_status <- function(value) {
  if (is.null(value)) return("skipped")
  if (inherits(value, "error")) return("failed")
  "success"
}

#' Return number of rows or missing numeric value
#' @keywords internal
nrow_or_na <- function(value) {
  if (is.null(value)) NA_real_ else nrow(value)
}

#' Normalize optional scalar to character
#' @keywords internal
scalar_character <- function(value) {
  if (is.null(value) || length(value) == 0L || is.na(value[1])) NA_character_ else as.character(value[1])
}

#' Normalize optional scalar to numeric
#' @keywords internal
scalar_numeric <- function(value) {
  if (is.null(value) || length(value) == 0L || is.na(value[1])) NA_real_ else as.numeric(value[1])
}

#' Return empty axis-model status
#' @keywords internal
empty_axis_model_status <- function() {
  data.frame(
    response_col = character(), status = character(), n_observations = numeric(),
    n_mesh_vertices = numeric(), waic = numeric(), dic = numeric(),
    error_message = character(), error_class = character(), stringsAsFactors = FALSE
  )
}

#' Return empty axis-prediction status
#' @keywords internal
empty_axis_prediction_status <- function() {
  data.frame(
    response_col = character(), status = character(), n_prediction_cells = numeric(),
    error_message = character(), error_class = character(), stringsAsFactors = FALSE
  )
}

#' Return empty workflow error log
#' @keywords internal
empty_workflow_error_log <- function() {
  data.frame(
    stage = character(), response_col = character(), message = character(),
    class = character(), stringsAsFactors = FALSE
  )
}
