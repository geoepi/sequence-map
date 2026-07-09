#' Convert predicted axis surfaces to a raster
#'
#' Converts projected point-centroid predictions into one or more terra raster
#' layers. These rasters are derived reporting products and do not change the
#' fitted SPDE model.
#'
#' @param surface Projected sf point predictions.
#' @param value_cols Optional prediction columns to rasterize.
#' @param resolution Optional positive raster cell size.
#' @param crs Optional output CRS override.
#' @param template Optional terra SpatRaster template.
#' @param fill Value for raster cells without a prediction.
#'
#' @return A terra SpatRaster.
#' @export
axis_surface_to_raster <- function(
  surface,
  value_cols = NULL,
  resolution = NULL,
  crs = NULL,
  template = NULL,
  fill = NA_real_
) {
  validate_axis_surface_points(surface)
  value_cols <- resolve_axis_surface_value_cols(surface, value_cols)

  if (!is.numeric(fill) || length(fill) != 1L) {
    stop("'fill' must be a single numeric value or NA_real_.", call. = FALSE)
  }

  raster_template <- prepare_axis_raster_template(
    surface = surface,
    resolution = resolution,
    crs = crs,
    template = template
  )
  surface <- sf::st_transform(surface, sf::st_crs(terra::crs(raster_template)))
  point_vector <- terra::vect(surface)
  raster_layers <- lapply(value_cols, function(value_col) {
    terra::rasterize(
      x = point_vector,
      y = raster_template,
      field = value_col,
      background = fill
    )
  })
  output_raster <- do.call(c, raster_layers)
  names(output_raster) <- value_cols
  output_raster
}

#' Write predicted axis surfaces as GeoTIFF files
#'
#' Writes one multiband GeoTIFF and one single-band GeoTIFF per selected axis
#' prediction or uncertainty column.
#'
#' @param surface An sf prediction surface or sequencemap_axis_surfaces object.
#' @param output_dir Directory for GeoTIFF output.
#' @param value_cols Optional prediction columns to rasterize.
#' @param prefix Output file-name prefix.
#' @param resolution Optional raster cell size.
#' @param template Optional terra SpatRaster template.
#' @param overwrite Logical; overwrite existing files.
#' @param datatype GDAL raster datatype.
#' @param gdal GDAL creation options.
#'
#' @return A list containing the raster and written file paths.
#' @export
write_axis_geotiffs <- function(
  surface,
  output_dir,
  value_cols = NULL,
  prefix = "sequence_map_axis_surface",
  resolution = NULL,
  template = NULL,
  overwrite = TRUE,
  datatype = "FLT4S",
  gdal = c("COMPRESS=LZW")
) {
  surface <- resolve_axis_surface(surface)
  if (!is.character(output_dir) || length(output_dir) != 1L || !nzchar(output_dir)) {
    stop("'output_dir' must be a non-empty directory path.", call. = FALSE)
  }
  if (!is.character(prefix) || length(prefix) != 1L || !nzchar(prefix)) {
    stop("'prefix' must be a non-empty character string.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  raster <- axis_surface_to_raster(
    surface = surface,
    value_cols = value_cols,
    resolution = resolution,
    template = template
  )
  value_cols <- names(raster)
  multiband_file <- file.path(output_dir, paste0(prefix, "_multiband.tif"))
  terra::writeRaster(
    raster, filename = multiband_file, overwrite = overwrite,
    datatype = datatype, gdal = gdal
  )
  singleband_files <- stats::setNames(
    vapply(value_cols, function(value_col) {
      filename <- file.path(output_dir, paste0(prefix, "_", value_col, ".tif"))
      terra::writeRaster(
        raster[[value_col]], filename = filename, overwrite = overwrite,
        datatype = datatype, gdal = gdal
      )
      filename
    }, character(1)),
    value_cols
  )
  list(
    raster = raster,
    multiband_file = multiband_file,
    singleband_files = singleband_files,
    value_cols = value_cols
  )
}

#' Aggregate axis rasters to polygons
#'
#' Calculates polygon-level reporting summaries from continuous predicted axis
#' rasters. These summaries are not separate polygon-level models and should
#' not be interpreted as nucleotide diversity.
#'
#' @param raster A terra SpatRaster or GeoTIFF path.
#' @param polygons An sf polygon object or path readable by sf.
#' @param polygon_id_col Optional polygon identifier column.
#' @param stats Statistics to calculate: "mean", "sd", "min", or "max".
#' @param na.rm Logical; remove missing raster values before summarizing.
#'
#' @return A data frame with one row per polygon.
#' @export
aggregate_axis_rasters_to_polygons <- function(
  raster,
  polygons,
  polygon_id_col = NULL,
  stats = c("mean", "sd", "min", "max"),
  na.rm = TRUE
) {
  raster <- resolve_axis_raster(raster)
  polygons <- resolve_axis_polygons(polygons)
  validate_polygon_summary_inputs(raster, polygons, polygon_id_col, stats, na.rm)

  polygons <- sf::st_transform(polygons, terra::crs(raster))
  extracted <- terra::extract(raster, terra::vect(polygons), ID = TRUE, exact = FALSE)
  layer_names <- names(raster)
  summaries <- data.frame(
    polygon_row = seq_len(nrow(polygons)),
    n_grid_cells = integer(nrow(polygons)),
    stringsAsFactors = FALSE
  )
  if (!is.null(polygon_id_col)) {
    summaries$polygon_id <- polygons[[polygon_id_col]]
  }

  for (polygon_index in seq_len(nrow(polygons))) {
    polygon_values <- extracted[extracted$ID == polygon_index, layer_names, drop = FALSE]
    summaries$n_grid_cells[polygon_index] <- nrow(polygon_values)
    for (layer_name in layer_names) {
      values <- polygon_values[[layer_name]]
      for (statistic in stats) {
        summaries[[paste0(layer_name, "_", statistic)]][polygon_index] <-
          calculate_polygon_raster_statistic(values, statistic, na.rm)
      }
    }
  }
  summaries
}

#' Validate a projected prediction surface
#'
#' @param surface An sf point surface.
#'
#' @return TRUE, invisibly.
#' @keywords internal
validate_axis_surface_points <- function(surface) {
  if (!inherits(surface, "sf")) {
    stop("'surface' must be an sf object containing point predictions.", call. = FALSE)
  }
  if (is.na(sf::st_crs(surface))) {
    stop("'surface' must have a valid projected CRS.", call. = FALSE)
  }
  if (isTRUE(sf::st_is_longlat(surface))) {
    stop("'surface' must use a projected CRS, not longitude/latitude.", call. = FALSE)
  }
  if (!all(as.character(sf::st_geometry_type(surface)) == "POINT")) {
    stop("'surface' must contain POINT geometries.", call. = FALSE)
  }
  if (nrow(surface) == 0L) {
    stop("'surface' must contain at least one point.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Resolve prediction columns for rasterization
#'
#' @param surface An sf surface.
#' @param value_cols Optional requested columns.
#'
#' @return Character vector.
#' @keywords internal
resolve_axis_surface_value_cols <- function(surface, value_cols = NULL) {
  if (is.null(value_cols)) {
    value_cols <- grep("(_pred|_sd|_q025|_q975)$", names(surface), value = TRUE)
  }
  if (!is.character(value_cols) || length(value_cols) == 0L) {
    stop("No prediction value columns were supplied or detected.", call. = FALSE)
  }
  missing_cols <- setdiff(value_cols, names(surface))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("'surface' is missing selected value columns: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  surface_data <- sf::st_drop_geometry(surface)
  nonnumeric_cols <- value_cols[!vapply(surface_data[value_cols], is.numeric, logical(1))]
  if (length(nonnumeric_cols) > 0L) {
    stop(
      sprintf("Selected raster value columns must be numeric: %s", paste(nonnumeric_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  value_cols
}

#' Prepare a raster template from point centroids
#'
#' @param surface An sf surface.
#' @param resolution Optional cell size.
#' @param crs Optional CRS override.
#' @param template Optional raster template.
#'
#' @return A terra SpatRaster.
#' @keywords internal
prepare_axis_raster_template <- function(surface, resolution = NULL, crs = NULL, template = NULL) {
  if (!is.null(template)) {
    if (!inherits(template, "SpatRaster")) {
      stop("'template' must be a terra SpatRaster.", call. = FALSE)
    }
    if (isTRUE(terra::is.lonlat(template))) {
      stop("'template' must use a projected CRS, not longitude/latitude.", call. = FALSE)
    }
    if (!is.null(crs) && sf::st_crs(crs) != sf::st_crs(terra::crs(template))) {
      stop("'crs' must match the CRS of 'template' when both are supplied.", call. = FALSE)
    }
    return(template[[1]])
  }

  output_crs <- if (is.null(crs)) sf::st_crs(surface) else sf::st_crs(crs)
  if (is.na(output_crs)) {
    stop("'crs' could not be resolved to a valid CRS.", call. = FALSE)
  }
  surface <- sf::st_transform(surface, output_crs)
  if (isTRUE(sf::st_is_longlat(surface))) {
    stop("Raster output must use a projected CRS, not longitude/latitude.", call. = FALSE)
  }

  if (is.null(resolution)) {
    coordinates <- sf::st_coordinates(surface)
    resolution <- c(
      infer_axis_grid_spacing(coordinates[, 1], "x"),
      infer_axis_grid_spacing(coordinates[, 2], "y")
    )
  } else {
    if (!is.numeric(resolution) || length(resolution) < 1L || length(resolution) > 2L ||
        any(is.na(resolution)) || any(resolution <= 0)) {
      stop("'resolution' must be a positive numeric value of length one or two.", call. = FALSE)
    }
    if (length(resolution) == 1L) {
      resolution <- rep(resolution, 2L)
    }
  }

  coordinates <- sf::st_coordinates(surface)
  terra::rast(
    xmin = min(coordinates[, 1]) - resolution[1] / 2,
    xmax = max(coordinates[, 1]) + resolution[1] / 2,
    ymin = min(coordinates[, 2]) - resolution[2] / 2,
    ymax = max(coordinates[, 2]) + resolution[2] / 2,
    resolution = resolution,
    crs = sf::st_crs(surface)$wkt
  )
}

#' Infer regular point-grid spacing
#'
#' @param values Numeric coordinate values.
#' @param dimension Coordinate label for error reporting.
#'
#' @return A positive numeric spacing.
#' @keywords internal
infer_axis_grid_spacing <- function(values, dimension) {
  spacing <- diff(sort(unique(values)))
  spacing <- spacing[spacing > 0]
  if (length(spacing) == 0L) {
    stop(
      sprintf("Could not infer %s resolution from point coordinates; supply 'resolution'.", dimension),
      call. = FALSE
    )
  }
  min(spacing)
}

#' Resolve a prediction surface from supported input forms
#'
#' @param surface Surface or batch surface object.
#'
#' @return An sf object.
#' @keywords internal
resolve_axis_surface <- function(surface) {
  if (inherits(surface, "sequencemap_axis_surfaces")) {
    surface <- surface$combined_surface
  }
  if (is.null(surface)) {
    stop("'surface' does not contain a combined prediction surface.", call. = FALSE)
  }
  surface
}

#' Resolve a raster input
#'
#' @param raster Raster or path.
#'
#' @return A terra SpatRaster.
#' @keywords internal
resolve_axis_raster <- function(raster) {
  if (inherits(raster, "SpatRaster")) {
    return(raster)
  }
  if (is.character(raster) && length(raster) == 1L && file.exists(raster)) {
    return(terra::rast(raster))
  }
  stop("'raster' must be a terra SpatRaster or an existing raster file path.", call. = FALSE)
}

#' Resolve polygon inputs
#'
#' @param polygons Sf polygons or a path.
#'
#' @return An sf object.
#' @keywords internal
resolve_axis_polygons <- function(polygons) {
  if (inherits(polygons, "sf")) {
    return(polygons)
  }
  if (is.character(polygons) && length(polygons) == 1L) {
    return(read_boundary(polygons))
  }
  stop("'polygons' must be an sf object or a readable boundary path.", call. = FALSE)
}

#' Validate polygon-raster summary inputs
#'
#' @param raster A raster.
#' @param polygons Polygon sf object.
#' @param polygon_id_col Optional identifier column.
#' @param stats Requested statistics.
#' @param na.rm Logical missing-value flag.
#'
#' @return TRUE, invisibly.
#' @keywords internal
validate_polygon_summary_inputs <- function(raster, polygons, polygon_id_col, stats, na.rm) {
  if (is.na(sf::st_crs(polygons))) {
    stop("'polygons' must have a valid CRS.", call. = FALSE)
  }
  polygon_types <- as.character(sf::st_geometry_type(polygons))
  if (!all(polygon_types %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("'polygons' must contain polygon geometries.", call. = FALSE)
  }
  if (!is.null(polygon_id_col) && (!is.character(polygon_id_col) || length(polygon_id_col) != 1L ||
      !polygon_id_col %in% names(polygons))) {
    stop("'polygon_id_col' must name a column in 'polygons'.", call. = FALSE)
  }
  valid_stats <- c("mean", "sd", "min", "max")
  if (!is.character(stats) || length(stats) == 0L || any(!stats %in% valid_stats)) {
    stop("'stats' must contain one or more of: mean, sd, min, max.", call. = FALSE)
  }
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    stop("'na.rm' must be a single non-missing logical value.", call. = FALSE)
  }
  if (!nzchar(terra::crs(raster))) {
    stop("'raster' must have a valid CRS.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Calculate one polygon-raster statistic
#'
#' @param values Raster cell values.
#' @param statistic Statistic name.
#' @param na.rm Logical missing-value flag.
#'
#' @return A numeric scalar.
#' @keywords internal
calculate_polygon_raster_statistic <- function(values, statistic, na.rm) {
  if (length(values) == 0L || (isTRUE(na.rm) && all(is.na(values)))) {
    return(NA_real_)
  }
  switch(
    statistic,
    mean = mean(values, na.rm = na.rm),
    sd = stats::sd(values, na.rm = na.rm),
    min = min(values, na.rm = na.rm),
    max = max(values, na.rm = na.rm)
  )
}
