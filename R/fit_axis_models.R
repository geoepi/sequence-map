#' Fit Gaussian SPDE models for multiple ordination axes
#'
#' Thin batch wrapper around [fit_axis_model()] for location-level ordination
#' responses such as `PC1_mean`, `PC2_mean`, `LD1_mean`, and `LD2_mean`.
#'
#' @param points_sf Projected `sf` points from [prepare_spatial_points()].
#' @param response_cols Optional character vector of response columns.
#' @param mesh Either an `INLA` mesh object or the list returned by
#'   [make_inla_mesh()].
#' @param iid_effects Optional character vector of iid random-effect columns.
#' @param family Likelihood family passed through to [fit_axis_model()].
#' @param compute_criteria Logical; passed through to [fit_axis_model()].
#' @param verbose Logical; passed through to [fit_axis_model()].
#' @param continue_on_error Logical; if `TRUE`, keep fitting after failures.
#'
#' @return A list of class `sequencemap_axis_models`.
#' @export
fit_axis_models <- function(
  points_sf,
  response_cols = NULL,
  mesh,
  iid_effects = NULL,
  family = "gaussian",
  compute_criteria = TRUE,
  verbose = FALSE,
  continue_on_error = TRUE
) {
  response_cols <- detect_axis_mean_response_cols(points_sf, response_cols = response_cols)

  models <- list()
  errors <- list()

  for (response_col in response_cols) {
    fit_result <- tryCatch(
      fit_axis_model(
        points_sf = points_sf,
        response_col = response_col,
        mesh = mesh,
        iid_effects = iid_effects,
        family = family,
        compute_criteria = compute_criteria,
        verbose = verbose
      ),
      error = function(error) error
    )

    if (inherits(fit_result, "error")) {
      if (!isTRUE(continue_on_error)) {
        stop(fit_result)
      }
      errors[[response_col]] <- make_batch_error_record(response_col, fit_result)
    } else {
      models[[response_col]] <- fit_result
    }
  }

  successful_response_cols <- names(models)
  failed_response_cols <- names(errors)

  model_summaries <- if (length(models) > 0L) {
    do.call(
      rbind,
      lapply(models, function(model_object) {
        data.frame(
          response_col = model_object$model_summary$response_col,
          n_observations = model_object$model_summary$n_observations,
          n_mesh_vertices = model_object$model_summary$n_mesh_vertices,
          waic = model_object$model_summary$waic,
          dic = model_object$model_summary$dic,
          stringsAsFactors = FALSE
        )
      })
    )
  } else {
    data.frame(
      response_col = character(0),
      n_observations = numeric(0),
      n_mesh_vertices = numeric(0),
      waic = numeric(0),
      dic = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      models = models,
      errors = errors,
      response_cols = response_cols,
      successful_response_cols = successful_response_cols,
      failed_response_cols = failed_response_cols,
      model_summaries = model_summaries,
      n_successful = length(successful_response_cols),
      n_failed = length(failed_response_cols)
    ),
    class = "sequencemap_axis_models"
  )
}

#' Detect batch response columns for axis-model fitting
#'
#' @param points_sf An sf object or data frame.
#' @param response_cols Optional explicit response columns.
#'
#' @return Character vector of response columns.
#' @keywords internal
detect_axis_mean_response_cols <- function(points_sf, response_cols = NULL) {
  if (!is.null(response_cols)) {
    missing_cols <- setdiff(response_cols, names(points_sf))
    if (length(missing_cols) > 0L) {
      stop(
        sprintf("Requested response columns were not found: %s", paste(missing_cols, collapse = ", ")),
        call. = FALSE
      )
    }
    return(response_cols)
  }

  candidate_cols <- grep("^(PC|PC_|LD|LD_).+_mean$", names(points_sf), value = TRUE)
  excluded_cols <- grep("(_sd|_min|_max|_q025|_q975|_pred)$|^(mean|sd|q025|q975)$", candidate_cols, value = TRUE)
  candidate_cols <- setdiff(candidate_cols, excluded_cols)
  candidate_cols <- candidate_cols[vapply(points_sf[candidate_cols], is.numeric, logical(1))]

  if (length(candidate_cols) == 0L) {
    stop("No valid ordination-axis mean response columns were found.", call. = FALSE)
  }

  candidate_cols
}

#' Convert an error to a structured batch record
#'
#' @param response_col Response column name.
#' @param error Error object.
#'
#' @return Named list.
#' @keywords internal
make_batch_error_record <- function(response_col, error) {
  list(
    response_col = response_col,
    message = conditionMessage(error),
    class = class(error)
  )
}
