#' Prepare projected spatial points for sequence mapping
#'
#' Converts a location summary table into a projected `sf` point object
#' suitable for downstream spatial modeling. Raw longitude/latitude input is
#' transformed to a projected CRS before returning.
#'
#' @param data A data frame containing coordinate columns.
#' @param lon_col Name of the x or longitude column.
#' @param lat_col Name of the y or latitude column.
#' @param crs_input Input CRS for the supplied coordinates.
#' @param crs_projected Target projected CRS. Required when `crs_input` is
#'   geographic.
#' @param id_col Optional identifier column to preserve.
#' @param keep_cols Optional character vector of additional columns to retain.
#' @param drop_missing Logical; if `TRUE`, drop rows with missing coordinates.
#'
#' @return A projected `sf` point object with numeric `x` and `y` columns.
#' @export
prepare_spatial_points <- function(
  data,
  lon_col = "longitude",
  lat_col = "latitude",
  crs_input = 4326,
  crs_projected = NULL,
  id_col = NULL,
  keep_cols = NULL,
  drop_missing = FALSE
) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  required_cols <- c(lon_col, lat_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("Coordinate columns were not found: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!is.null(id_col) && !id_col %in% names(data)) {
    stop(sprintf("Identifier column '%s' was not found.", id_col), call. = FALSE)
  }

  keep_cols <- unique(c(id_col, keep_cols, lon_col, lat_col))
  keep_cols <- keep_cols[!is.na(keep_cols)]
  keep_cols <- intersect(keep_cols, names(data))
  data_subset <- if (length(keep_cols) > 0L) {
    data[, keep_cols, drop = FALSE]
  } else {
    data
  }

  data_subset[[lon_col]] <- coerce_numeric_column(data_subset[[lon_col]], lon_col)
  data_subset[[lat_col]] <- coerce_numeric_column(data_subset[[lat_col]], lat_col)

  missing_coordinate_rows <- is.na(data_subset[[lon_col]]) | is.na(data_subset[[lat_col]])
  if (any(missing_coordinate_rows)) {
    if (!isTRUE(drop_missing)) {
      stop(
        "Missing coordinates detected. Set `drop_missing = TRUE` to discard those rows.",
        call. = FALSE
      )
    }
    data_subset <- data_subset[!missing_coordinate_rows, , drop = FALSE]
  }

  if (nrow(data_subset) == 0L) {
    stop("No rows remain after coordinate filtering.", call. = FALSE)
  }

  input_crs <- sf::st_crs(crs_input)
  if (is.na(input_crs)) {
    stop("`crs_input` could not be resolved to a valid CRS.", call. = FALSE)
  }

  points_sf <- sf::st_as_sf(
    data_subset,
    coords = c(lon_col, lat_col),
    crs = input_crs,
    remove = FALSE
  )

  input_is_longlat <- isTRUE(sf::st_is_longlat(points_sf))
  if (input_is_longlat && is.null(crs_projected)) {
    stop(
      "`crs_projected` is required when input coordinates are longitude/latitude.",
      call. = FALSE
    )
  }

  if (!is.null(crs_projected)) {
    projected_crs <- sf::st_crs(crs_projected)
    if (is.na(projected_crs)) {
      stop("`crs_projected` could not be resolved to a valid CRS.", call. = FALSE)
    }
    points_sf <- sf::st_transform(points_sf, projected_crs)
  } else {
    projected_crs <- sf::st_crs(points_sf)
  }

  if (isTRUE(sf::st_is_longlat(points_sf))) {
    stop("Spatial points must be projected before analysis.", call. = FALSE)
  }

  coords <- sf::st_coordinates(points_sf)
  points_sf$x <- coords[, 1]
  points_sf$y <- coords[, 2]
  points_sf
}
