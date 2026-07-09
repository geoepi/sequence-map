#' Create a prediction grid inside a boundary
#'
#' Builds a regular grid of point centroids inside a supplied boundary using a
#' projected CRS.
#'
#' @param boundary An `sf` object or a path readable by [sf::st_read()].
#' @param crs_projected Projected CRS used for analysis.
#' @param grid_resolution Numeric grid spacing in projected CRS units.
#' @param include_boundary Logical; if `TRUE`, keep centroids intersecting the
#'   boundary. If `FALSE`, require points to fall strictly within the polygon.
#' @param return Output type, either `"sf"` or `"terra"`.
#'
#' @return An `sf` object or `terra` vector of prediction points.
#' @export
make_prediction_grid <- function(
  boundary,
  crs_projected,
  grid_resolution,
  include_boundary = TRUE,
  return = c("sf", "terra")
) {
  return <- match.arg(return)
  boundary_sf <- normalize_boundary(boundary)

  if (missing(grid_resolution) || length(grid_resolution) != 1L || !is.numeric(grid_resolution) ||
      is.na(grid_resolution) || grid_resolution <= 0) {
    stop("`grid_resolution` must be a single positive numeric value.", call. = FALSE)
  }

  projected_crs <- sf::st_crs(crs_projected)
  if (is.na(projected_crs)) {
    stop("`crs_projected` could not be resolved to a valid CRS.", call. = FALSE)
  }

  boundary_sf <- sf::st_make_valid(boundary_sf)
  boundary_sf <- sf::st_transform(boundary_sf, projected_crs)

  if (nrow(boundary_sf) == 0L || all(sf::st_is_empty(boundary_sf))) {
    stop("`boundary` must contain non-empty geometries.", call. = FALSE)
  }

  bbox <- sf::st_bbox(boundary_sf)
  width <- as.numeric(bbox["xmax"] - bbox["xmin"])
  height <- as.numeric(bbox["ymax"] - bbox["ymin"])
  if (width <= 0 || height <= 0) {
    stop("`boundary` must span a positive area in projected coordinates.", call. = FALSE)
  }

  estimated_points <- ceiling(width / grid_resolution) * ceiling(height / grid_resolution)
  if (estimated_points > 1e6) {
    stop(
      "`grid_resolution` is too small for this boundary and would generate an excessively large grid.",
      call. = FALSE
    )
  }
  if (grid_resolution > max(width, height) * 100) {
    stop(
      "`grid_resolution` is too large relative to the boundary extent to generate useful grid cells.",
      call. = FALSE
    )
  }

  grid_geometry <- sf::st_make_grid(
    boundary_sf,
    cellsize = grid_resolution,
    what = "centers"
  )

  if (length(grid_geometry) == 0L) {
    stop("No grid points were generated from the requested boundary and resolution.", call. = FALSE)
  }

  grid_sf <- sf::st_sf(geometry = grid_geometry)
  inside <- if (isTRUE(include_boundary)) {
    lengths(sf::st_intersects(grid_sf, boundary_sf)) > 0L
  } else {
    lengths(sf::st_within(grid_sf, boundary_sf)) > 0L
  }
  grid_sf <- grid_sf[inside, , drop = FALSE]

  if (nrow(grid_sf) == 0L) {
    stop("No grid points fall inside the supplied boundary.", call. = FALSE)
  }

  coords <- sf::st_coordinates(grid_sf)
  grid_sf$grid_id <- sprintf("grid_%s", seq_len(nrow(grid_sf)))
  grid_sf$x <- coords[, 1]
  grid_sf$y <- coords[, 2]
  grid_sf <- grid_sf[, c("grid_id", "x", "y", "geometry")]
  sf::st_crs(grid_sf) <- projected_crs

  if (return == "terra") {
    return(terra::vect(grid_sf))
  }

  grid_sf
}
