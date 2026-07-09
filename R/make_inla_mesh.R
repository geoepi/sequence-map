#' Create an INLA mesh for sequence mapping
#'
#' Builds an `INLA` mesh from projected sampled points, optionally using a
#' boundary and prediction grid to define the spatial domain extent.
#'
#' @param points_sf Projected `sf` points returned by [prepare_spatial_points()].
#' @param boundary Optional boundary as an `sf` object or path readable by
#'   [sf::st_read()].
#' @param prediction_grid Optional projected prediction grid.
#' @param max_edge Numeric scalar or length-2 vector passed to
#'   `INLA::inla.mesh.2d()`.
#' @param cutoff Numeric cutoff passed to `INLA::inla.mesh.2d()`.
#' @param offset Numeric scalar or length-2 vector passed to
#'   `INLA::inla.mesh.2d()`.
#' @param crs_projected Optional projected CRS used to validate inputs.
#'
#' @return A list with mesh object, coordinates, CRS, point count, and a mesh
#'   summary table.
#' @export
make_inla_mesh <- function(
  points_sf,
  boundary = NULL,
  prediction_grid = NULL,
  max_edge,
  cutoff,
  offset,
  crs_projected = NULL
) {
  validate_projected_sf_points(points_sf, object_name = "points_sf")

  if (!missing(crs_projected) && !is.null(crs_projected)) {
    target_crs <- sf::st_crs(crs_projected)
    if (is.na(target_crs)) {
      stop("`crs_projected` could not be resolved to a valid CRS.", call. = FALSE)
    }
    if (sf::st_crs(points_sf) != target_crs) {
      stop("`points_sf` CRS does not match `crs_projected`.", call. = FALSE)
    }
  }

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required for `make_inla_mesh()`.", call. = FALSE)
  }

  if (missing(max_edge) || !is.numeric(max_edge) || length(max_edge) < 1L || length(max_edge) > 2L ||
      any(is.na(max_edge)) || any(max_edge <= 0)) {
    stop("`max_edge` must be a positive numeric vector of length 1 or 2.", call. = FALSE)
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff) || cutoff < 0) {
    stop("`cutoff` must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.numeric(offset) || length(offset) < 1L || length(offset) > 2L ||
      any(is.na(offset)) || any(offset <= 0)) {
    stop("`offset` must be a positive numeric vector of length 1 or 2.", call. = FALSE)
  }

  loc <- sf::st_coordinates(points_sf)[, 1:2, drop = FALSE]
  domain_coords <- loc
  crs <- sf::st_crs(points_sf)

  if (!is.null(prediction_grid)) {
    validate_projected_sf_points(prediction_grid, object_name = "prediction_grid")
    if (sf::st_crs(prediction_grid) != crs) {
      prediction_grid <- sf::st_transform(prediction_grid, crs)
    }
    domain_coords <- rbind(domain_coords, sf::st_coordinates(prediction_grid)[, 1:2, drop = FALSE])
  }

  boundary_segment <- NULL
  if (!is.null(boundary)) {
    boundary_sf <- normalize_boundary(boundary)
    boundary_sf <- sf::st_transform(sf::st_make_valid(boundary_sf), crs)
    boundary_coords <- extract_boundary_coordinates(boundary_sf)
    if (nrow(boundary_coords) < 3L) {
      stop("`boundary` must contain polygonal coordinates suitable for mesh construction.", call. = FALSE)
    }
    domain_coords <- rbind(domain_coords, boundary_coords)
    boundary_segment <- INLA::inla.nonconvex.hull(boundary_coords, convex = max(max_edge) * 2)
  } else if (nrow(domain_coords) >= 3L) {
    boundary_segment <- INLA::inla.nonconvex.hull(domain_coords, convex = max(max_edge) * 2)
  }

  mesh <- INLA::inla.mesh.2d(
    loc = loc,
    boundary = boundary_segment,
    max.edge = max_edge,
    cutoff = cutoff,
    offset = offset
  )

  vertex_loc <- mesh$loc[, 1:2, drop = FALSE]
  mesh_summary <- data.frame(
    n_vertices = nrow(vertex_loc),
    n_triangles = nrow(mesh$graph$tv),
    max_edge = paste(max_edge, collapse = ","),
    cutoff = cutoff,
    x_min = min(vertex_loc[, 1]),
    x_max = max(vertex_loc[, 1]),
    y_min = min(vertex_loc[, 2]),
    y_max = max(vertex_loc[, 2]),
    stringsAsFactors = FALSE
  )

  list(
    mesh = mesh,
    loc = loc,
    crs = crs,
    n_points = nrow(loc),
    mesh_summary = mesh_summary
  )
}

#' Normalize a boundary input to sf
#'
#' @param boundary A boundary object or path.
#'
#' @return An sf object.
#' @keywords internal
normalize_boundary <- function(boundary) {
  if (inherits(boundary, "sf")) {
    return(boundary)
  }

  if (is.character(boundary) && length(boundary) == 1L) {
    return(read_boundary(boundary))
  }

  stop("`boundary` must be an sf object or a readable boundary path.", call. = FALSE)
}

#' Validate projected sf point input
#'
#' @param x An sf object.
#' @param object_name Name used in error messages.
#'
#' @return `TRUE`, invisibly.
#' @keywords internal
validate_projected_sf_points <- function(x, object_name = "x") {
  if (!inherits(x, "sf")) {
    stop(sprintf("`%s` must be an sf object.", object_name), call. = FALSE)
  }

  if (is.na(sf::st_crs(x))) {
    stop(sprintf("`%s` must have a valid CRS.", object_name), call. = FALSE)
  }

  if (isTRUE(sf::st_is_longlat(x))) {
    stop(
      sprintf("`%s` appears to be longitude/latitude. Transform to a projected CRS first.", object_name),
      call. = FALSE
    )
  }

  geometry_types <- unique(as.character(sf::st_geometry_type(x)))
  if (!all(geometry_types %in% c("POINT", "MULTIPOINT"))) {
    stop(sprintf("`%s` must contain point geometries.", object_name), call. = FALSE)
  }

  invisible(TRUE)
}

#' Extract polygon boundary coordinates for mesh construction
#'
#' @param boundary_sf An sf boundary object.
#'
#' @return A two-column numeric matrix.
#' @keywords internal
extract_boundary_coordinates <- function(boundary_sf) {
  boundary_geom <- sf::st_union(boundary_sf)
  boundary_line <- sf::st_boundary(boundary_geom)
  boundary_points <- sf::st_cast(boundary_line, "POINT", warn = FALSE)
  coords <- sf::st_coordinates(boundary_points)

  if (nrow(coords) == 0L) {
    return(matrix(numeric(0), ncol = 2))
  }

  coords[, 1:2, drop = FALSE]
}
