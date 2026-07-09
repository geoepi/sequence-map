test_that("prepare_spatial_points transforms lon lat to projected coordinates", {
  data <- data.frame(
    sequence_id = c("a", "b"),
    longitude = c(-93.1, -93.2),
    latitude = c(44.9, 45.0),
    value = c(1, 2),
    stringsAsFactors = FALSE
  )

  result <- prepare_spatial_points(
    data,
    crs_input = 4326,
    crs_projected = 3857,
    id_col = "sequence_id",
    keep_cols = "value"
  )

  expect_true(inherits(result, "sf"))
  expect_true(all(c("x", "y", "sequence_id", "value") %in% names(result)))
  expect_false(any(abs(result$x - data$longitude) < 1e-6))
  expect_false(isTRUE(sf::st_is_longlat(result)))
})

test_that("prepare_spatial_points errors on missing coordinates unless drop_missing is true", {
  data <- data.frame(
    longitude = c(-93.1, NA_real_),
    latitude = c(44.9, 45.0),
    stringsAsFactors = FALSE
  )

  expect_error(
    prepare_spatial_points(data, crs_input = 4326, crs_projected = 3857),
    "Missing coordinates detected"
  )

  result <- prepare_spatial_points(
    data,
    crs_input = 4326,
    crs_projected = 3857,
    drop_missing = TRUE
  )

  expect_equal(nrow(result), 1)
})

test_that("make_prediction_grid creates grid points inside a simple polygon", {
  boundary <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)
      ))),
      crs = 3857
    )
  )

  grid <- make_prediction_grid(
    boundary = boundary,
    crs_projected = 3857,
    grid_resolution = 500
  )

  expect_true(inherits(grid, "sf"))
  expect_true(nrow(grid) > 0)
  expect_true(all(c("grid_id", "x", "y") %in% names(grid)))
})

test_that("make_prediction_grid returns grid identifiers coordinates and geometry", {
  boundary <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)
      ))),
      crs = 3857
    )
  )

  grid <- make_prediction_grid(boundary, crs_projected = 3857, grid_resolution = 400)

  expect_true(all(c("grid_id", "x", "y", "geometry") %in% names(grid)))
})

test_that("make_prediction_grid errors when no grid points are produced", {
  boundary <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(0, 0), c(1000, 0), c(1000, 100), c(0, 100), c(0, 0)
      ))),
      crs = 3857
    )
  )

  expect_error(
    make_prediction_grid(
      boundary,
      crs_projected = 3857,
      grid_resolution = 500,
      include_boundary = FALSE
    ),
    "No grid points fall inside"
  )
})

test_that("make_inla_mesh skips gracefully if INLA is not installed", {
  if (requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is installed; skip not-installed check.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(x = c(0, 1), y = c(0, 1)),
    coords = c("x", "y"),
    crs = 3857
  )

  expect_error(
    make_inla_mesh(points_sf, max_edge = c(100, 200), cutoff = 10, offset = c(50, 100)),
    "Package 'INLA' is required"
  )
})

test_that("make_inla_mesh returns mesh components when INLA is available", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(id = 1:4, x = c(0, 100, 0, 100), y = c(0, 0, 100, 100)),
    coords = c("x", "y"),
    crs = 3857
  )
  boundary <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(-50, -50), c(150, -50), c(150, 150), c(-50, 150), c(-50, -50)
      ))),
      crs = 3857
    )
  )

  mesh_result <- make_inla_mesh(
    points_sf = points_sf,
    boundary = boundary,
    max_edge = c(50, 100),
    cutoff = 1,
    offset = c(20, 40),
    crs_projected = 3857
  )

  expect_true(all(c("mesh", "loc", "crs", "n_points", "mesh_summary") %in% names(mesh_result)))
  expect_true(is.matrix(mesh_result$loc))
  expect_equal(mesh_result$n_points, 4)
})

test_that("make_inla_mesh errors on unprojected longitude latitude points", {
  points_sf <- sf::st_as_sf(
    data.frame(longitude = c(-93.1, -93.0), latitude = c(44.9, 45.0)),
    coords = c("longitude", "latitude"),
    crs = 4326
  )

  expect_error(
    make_inla_mesh(points_sf, max_edge = c(100, 200), cutoff = 10, offset = c(50, 100)),
    "longitude/latitude"
  )
})
