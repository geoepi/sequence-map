make_test_axis_surface <- function() {
  grid <- expand.grid(x = c(0, 100), y = c(0, 100))
  sf::st_as_sf(
    transform(
      grid,
      PC1_mean_pred = c(1, 2, 3, 4),
      PC1_mean_sd = c(0.1, 0.2, 0.3, 0.4),
      label = c("a", "b", "c", "d")
    ),
    coords = c("x", "y"),
    crs = 3857,
    remove = FALSE
  )
}

make_test_polygons <- function() {
  sf::st_sf(
    province = c("west", "east"),
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(-50, -50), c(50, -50), c(50, 150), c(-50, 150), c(-50, -50)))),
      sf::st_polygon(list(rbind(c(50, -50), c(150, -50), c(150, 150), c(50, 150), c(50, -50)))),
      crs = 3857
    )
  )
}

test_that("axis_surface_to_raster converts a regular point grid", {
  raster <- axis_surface_to_raster(make_test_axis_surface(), value_cols = "PC1_mean_pred")

  expect_s4_class(raster, "SpatRaster")
  expect_equal(terra::nlyr(raster), 1)
  expect_equal(names(raster), "PC1_mean_pred")
  expect_equal(terra::res(raster), c(100, 100))
})

test_that("axis_surface_to_raster detects prediction and uncertainty columns", {
  raster <- axis_surface_to_raster(make_test_axis_surface())

  expect_equal(
    names(raster),
    c("PC1_mean_pred", "PC1_mean_sd")
  )
})

test_that("axis_surface_to_raster rejects unprojected input", {
  surface <- sf::st_transform(make_test_axis_surface(), 4326)

  expect_error(
    axis_surface_to_raster(surface),
    "projected CRS"
  )
})

test_that("axis_surface_to_raster rejects nonnumeric selected values", {
  expect_error(
    axis_surface_to_raster(make_test_axis_surface(), value_cols = "label"),
    "must be numeric"
  )
})

test_that("write_axis_geotiffs writes multiband and single-band files", {
  output_dir <- file.path(tempdir(), "sequencemap-geotiffs")
  result <- write_axis_geotiffs(make_test_axis_surface(), output_dir = output_dir)

  expect_s4_class(result$raster, "SpatRaster")
  expect_true(file.exists(result$multiband_file))
  expect_true(all(file.exists(result$singleband_files)))
  expect_equal(length(result$singleband_files), 2)
})

test_that("write_axis_geotiffs accepts a batch prediction object", {
  surface_object <- structure(
    list(combined_surface = make_test_axis_surface()),
    class = "sequencemap_axis_surfaces"
  )
  result <- write_axis_geotiffs(
    surface_object,
    output_dir = file.path(tempdir(), "sequencemap-geotiffs-batch")
  )

  expect_true(file.exists(result$multiband_file))
})

test_that("aggregate_axis_rasters_to_polygons returns polygon means", {
  raster <- axis_surface_to_raster(make_test_axis_surface())
  summaries <- aggregate_axis_rasters_to_polygons(
    raster,
    make_test_polygons(),
    stats = "mean"
  )

  expect_equal(nrow(summaries), 2)
  expect_true("PC1_mean_pred_mean" %in% names(summaries))
  expect_equal(summaries$PC1_mean_pred_mean, c(2, 3))
  expect_true(all(summaries$n_grid_cells == 2))
})

test_that("aggregate_axis_rasters_to_polygons preserves polygon identifiers", {
  summaries <- aggregate_axis_rasters_to_polygons(
    axis_surface_to_raster(make_test_axis_surface()),
    make_test_polygons(),
    polygon_id_col = "province",
    stats = "mean"
  )

  expect_equal(summaries$polygon_id, c("west", "east"))
})

test_that("aggregate_axis_rasters_to_polygons transforms polygon CRS", {
  polygons_longlat <- sf::st_transform(make_test_polygons(), 4326)
  summaries <- aggregate_axis_rasters_to_polygons(
    axis_surface_to_raster(make_test_axis_surface()),
    polygons_longlat,
    polygon_id_col = "province",
    stats = "mean"
  )

  expect_equal(summaries$PC1_mean_pred_mean, c(2, 3))
})
