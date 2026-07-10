#' Check spatial workflow dependencies
#'
#' @return TRUE, invisibly.
#' @keywords internal
check_spatial_workflow_dependencies <- function() {
  required_packages <- c("sf", "terra", "INLA")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop(
      sprintf(
        "Spatial modeling cannot start because required package(s) are unavailable: %s. Install missing packages and retry.",
        paste(missing_packages, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Check location support before spatial modeling
#'
#' @param points_sf Projected sampled locations.
#'
#' @return TRUE, invisibly.
#' @keywords internal
check_spatial_location_support <- function(points_sf) {
  if (nrow(points_sf) < 3L) {
    stop(
      sprintf(
        "Spatial modeling requires at least 3 sampled locations after aggregation; found %d.",
        nrow(points_sf)
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
