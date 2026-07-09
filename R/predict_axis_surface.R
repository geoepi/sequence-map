#' Predict an ordination axis surface on a spatial grid
#'
#' Projects a fitted single-axis SPDE model to a prediction grid and returns
#' posterior mean and approximate uncertainty for each grid cell.
#'
#' @param axis_model Object returned by [fit_axis_model()].
#' @param prediction_grid Projected `sf` grid from [make_prediction_grid()].
#' @param include_observed Logical; reserved for future use.
#'
#' @return The prediction grid with added posterior summary columns.
#' @export
predict_axis_surface <- function(
  axis_model,
  prediction_grid,
  include_observed = FALSE
) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required for `predict_axis_surface()`.", call. = FALSE)
  }

  required_components <- c("model", "spde", "mesh", "response_col", "fitted_data")
  missing_components <- setdiff(required_components, names(axis_model))
  if (length(missing_components) > 0L) {
    stop(
      sprintf(
        "`axis_model` is missing required components: %s",
        paste(missing_components, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  validate_projected_sf_points(prediction_grid, object_name = "prediction_grid")

  model_crs <- sf::st_crs(axis_model$fitted_data)
  if (sf::st_crs(prediction_grid) != model_crs) {
    stop("`prediction_grid` CRS must match the fitted model CRS.", call. = FALSE)
  }

  response_col <- axis_model$response_col
  grid_coords <- sf::st_coordinates(prediction_grid)[, 1:2, drop = FALSE]
  a_pred <- INLA::inla.spde.make.A(mesh = axis_model$mesh, loc = grid_coords)

  intercept_summary <- axis_model$model$summary.fixed
  if (!"intercept" %in% rownames(intercept_summary)) {
    stop("Fitted model is missing the intercept summary required for prediction.", call. = FALSE)
  }

  spatial_summary <- axis_model$model$summary.random$spatial_field
  if (is.null(spatial_summary)) {
    stop("Fitted model is missing the spatial field summary required for prediction.", call. = FALSE)
  }

  intercept_mean <- intercept_summary["intercept", "mean"]
  intercept_sd <- intercept_summary["intercept", "sd"]
  spatial_mean <- spatial_summary$mean
  spatial_sd <- spatial_summary$sd

  mean_values <- as.numeric(intercept_mean + as.vector(a_pred %*% spatial_mean))
  spatial_var <- as.numeric((a_pred^2) %*% (spatial_sd^2))
  sd_values <- sqrt(intercept_sd^2 + pmax(spatial_var, 0))
  q025_values <- mean_values - 1.96 * sd_values
  q975_values <- mean_values + 1.96 * sd_values

  prediction_sf <- prediction_grid
  prediction_sf$response_col <- response_col
  prediction_sf$mean <- mean_values
  prediction_sf$sd <- sd_values
  prediction_sf$q025 <- q025_values
  prediction_sf$q975 <- q975_values

  prediction_sf[[paste0(response_col, "_pred")]] <- mean_values
  prediction_sf[[paste0(response_col, "_sd")]] <- sd_values
  prediction_sf[[paste0(response_col, "_q025")]] <- q025_values
  prediction_sf[[paste0(response_col, "_q975")]] <- q975_values

  if (isTRUE(include_observed)) {
    prediction_sf$include_observed <- TRUE
  }

  prediction_sf
}
