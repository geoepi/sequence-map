#' Predict multiple axis surfaces on a common grid
#'
#' Thin batch wrapper around [predict_axis_surface()] for a set of fitted
#' location-level axis models.
#'
#' @param axis_models Object returned by [fit_axis_models()] or a named list of
#'   [fit_axis_model()] results.
#' @param prediction_grid Projected `sf` grid from [make_prediction_grid()].
#' @param include_observed Logical; passed through to [predict_axis_surface()].
#' @param continue_on_error Logical; if `TRUE`, keep predicting after failures.
#'
#' @return A list of class `sequencemap_axis_surfaces`.
#' @export
predict_axis_surfaces <- function(
  axis_models,
  prediction_grid,
  include_observed = FALSE,
  continue_on_error = TRUE
) {
  model_list <- normalize_axis_model_list(axis_models)
  response_cols <- names(model_list)

  surfaces <- list()
  errors <- list()

  for (response_col in response_cols) {
    prediction_result <- tryCatch(
      predict_axis_surface(
        axis_model = model_list[[response_col]],
        prediction_grid = prediction_grid,
        include_observed = include_observed
      ),
      error = function(error) error
    )

    if (inherits(prediction_result, "error")) {
      if (!isTRUE(continue_on_error)) {
        stop(prediction_result)
      }
      errors[[response_col]] <- make_batch_error_record(response_col, prediction_result)
    } else {
      surfaces[[response_col]] <- prediction_result
    }
  }

  successful_response_cols <- names(surfaces)
  failed_response_cols <- names(errors)

  combined_surface <- combine_axis_surfaces(surfaces)

  structure(
    list(
      surfaces = surfaces,
      combined_surface = combined_surface,
      errors = errors,
      response_cols = response_cols,
      successful_response_cols = successful_response_cols,
      failed_response_cols = failed_response_cols
    ),
    class = "sequencemap_axis_surfaces"
  )
}

#' Normalize a batch axis-model input to a named list
#'
#' @param axis_models Batch object or named list.
#'
#' @return Named list of model objects.
#' @keywords internal
normalize_axis_model_list <- function(axis_models) {
  if (inherits(axis_models, "sequencemap_axis_models")) {
    return(axis_models$models)
  }

  if (!is.list(axis_models) || length(axis_models) == 0L) {
    stop("`axis_models` must be a non-empty batch object or named list of axis models.", call. = FALSE)
  }

  model_names <- names(axis_models)
  if (is.null(model_names) || anyNA(model_names) || any(!nzchar(model_names))) {
    model_names <- vapply(axis_models, function(model_object) {
      if (!is.null(model_object$response_col)) model_object$response_col else ""
    }, character(1))
    if (any(!nzchar(model_names))) {
      stop("`axis_models` must be named, or each model must contain `response_col`.", call. = FALSE)
    }
    names(axis_models) <- model_names
  }

  axis_models
}

#' Combine response-specific surface outputs onto one grid
#'
#' @param surfaces Named list of sf surfaces.
#'
#' @return An sf object or `NULL`.
#' @keywords internal
combine_axis_surfaces <- function(surfaces) {
  if (length(surfaces) == 0L) {
    return(NULL)
  }

  first_surface <- surfaces[[1]]
  combined_surface <- first_surface[, intersect(c("grid_id", "x", "y", "geometry"), names(first_surface))]

  for (response_col in names(surfaces)) {
    surface <- surfaces[[response_col]]
    specific_cols <- intersect(
      c(
        paste0(response_col, "_pred"),
        paste0(response_col, "_sd"),
        paste0(response_col, "_q025"),
        paste0(response_col, "_q975")
      ),
      names(surface)
    )

    if (length(specific_cols) == 0L) {
      next
    }

    surface_data <- sf::st_drop_geometry(surface)[, c("grid_id", specific_cols), drop = FALSE]
    combined_surface <- merge(
      combined_surface,
      surface_data,
      by = "grid_id",
      all.x = TRUE,
      sort = FALSE
    )
    combined_surface <- sf::st_as_sf(combined_surface)
    sf::st_geometry(combined_surface) <- sf::st_geometry(first_surface)[match(combined_surface$grid_id, first_surface$grid_id)]
    sf::st_crs(combined_surface) <- sf::st_crs(first_surface)
  }

  combined_surface
}
