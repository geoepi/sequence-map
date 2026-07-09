test_that("axis model helpers skip gracefully if INLA is not installed", {
  if (requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is installed; skip not-installed branch.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(PC1_mean = c(0, 1, 2), x = c(0, 100, 200), y = c(0, 100, 0)),
    coords = c("x", "y"),
    crs = 3857
  )

  expect_error(
    fit_axis_model(
      points_sf = points_sf,
      response_col = "PC1_mean",
      mesh = list(mesh = list(loc = matrix(0, nrow = 3, ncol = 2)))
    ),
    "Package 'INLA' is required"
  )
})

test_that("fit_axis_model fits a single-axis projected model", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(
      location_id = paste0("loc", 1:4),
      PC1_mean = c(0.1, 0.3, 0.8, 1.0),
      x = c(0, 1000, 0, 1000),
      y = c(0, 0, 1000, 1000)
    ),
    coords = c("x", "y"),
    crs = 3857
  )
  mesh_obj <- make_inla_mesh(
    points_sf = points_sf,
    max_edge = c(800, 1600),
    cutoff = 1,
    offset = c(200, 400),
    crs_projected = 3857
  )

  result <- fit_axis_model(
    points_sf = points_sf,
    response_col = "PC1_mean",
    mesh = mesh_obj
  )

  expect_true(all(c("model", "spde", "mesh", "stack", "response_col", "fitted_data", "formula", "model_summary") %in% names(result)))
  expect_equal(result$response_col, "PC1_mean")
  expect_equal(result$model_summary$response_col, "PC1_mean")
  expect_equal(result$model_summary$n_observations, 4)
})

test_that("fit_axis_model errors for missing response column", {
  points_sf <- sf::st_as_sf(
    data.frame(x = c(0, 1000, 0), y = c(0, 0, 1000)),
    coords = c("x", "y"),
    crs = 3857
  )

  expect_error(
    fit_axis_model(points_sf, "PC1_mean", mesh = list(mesh = list(loc = matrix(0, nrow = 3, ncol = 2)))),
    "Response column 'PC1_mean' was not found"
  )
})

test_that("fit_axis_model errors for nonnumeric response", {
  points_sf <- sf::st_as_sf(
    data.frame(PC1_mean = c("a", "b", "c"), x = c(0, 1000, 0), y = c(0, 0, 1000)),
    coords = c("x", "y"),
    crs = 3857
  )

  expect_error(
    fit_axis_model(points_sf, "PC1_mean", mesh = list(mesh = list(loc = matrix(0, nrow = 3, ncol = 2)))),
    "must be numeric"
  )
})

test_that("fit_axis_model errors when fewer than three non-missing observations remain", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(PC1_mean = c(0.1, NA, NA), x = c(0, 1000, 0), y = c(0, 0, 1000)),
    coords = c("x", "y"),
    crs = 3857
  )
  mesh_obj <- make_inla_mesh(
    points_sf = sf::st_as_sf(
      data.frame(PC1_mean = c(0.1, 0.2, 0.3), x = c(0, 1000, 0), y = c(0, 0, 1000)),
      coords = c("x", "y"),
      crs = 3857
    ),
    max_edge = c(800, 1600),
    cutoff = 1,
    offset = c(200, 400),
    crs_projected = 3857
  )

  expect_warning(
    expect_error(
      fit_axis_model(points_sf, "PC1_mean", mesh = mesh_obj),
      "At least 3 non-missing observations"
    ),
    "Dropped 2 observations"
  )
})

test_that("predict_axis_surface returns posterior summaries on a grid", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(
      PC1_mean = c(0.1, 0.3, 0.8, 1.0),
      x = c(0, 1000, 0, 1000),
      y = c(0, 0, 1000, 1000)
    ),
    coords = c("x", "y"),
    crs = 3857
  )
  boundary <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(
        c(-200, -200), c(1200, -200), c(1200, 1200), c(-200, 1200), c(-200, -200)
      ))),
      crs = 3857
    )
  )
  grid <- make_prediction_grid(boundary, crs_projected = 3857, grid_resolution = 700)
  mesh_obj <- make_inla_mesh(
    points_sf = points_sf,
    boundary = boundary,
    prediction_grid = grid,
    max_edge = c(800, 1600),
    cutoff = 1,
    offset = c(200, 400),
    crs_projected = 3857
  )
  model_fit <- fit_axis_model(points_sf, "PC1_mean", mesh = mesh_obj)

  surface <- predict_axis_surface(model_fit, grid)

  expect_true(all(c("mean", "sd", "q025", "q975", "PC1_mean_pred", "PC1_mean_sd", "PC1_mean_q025", "PC1_mean_q975") %in% names(surface)))
  expect_equal(nrow(surface), nrow(grid))
})

test_that("predict_axis_surface errors on unprojected longitude latitude grid", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(
      PC1_mean = c(0.1, 0.3, 0.8, 1.0),
      x = c(0, 1000, 0, 1000),
      y = c(0, 0, 1000, 1000)
    ),
    coords = c("x", "y"),
    crs = 3857
  )
  mesh_obj <- make_inla_mesh(
    points_sf = points_sf,
    max_edge = c(800, 1600),
    cutoff = 1,
    offset = c(200, 400),
    crs_projected = 3857
  )
  model_fit <- fit_axis_model(points_sf, "PC1_mean", mesh = mesh_obj)
  grid_lonlat <- sf::st_as_sf(
    data.frame(longitude = c(-93.1, -93.0), latitude = c(44.9, 45.0)),
    coords = c("longitude", "latitude"),
    crs = 4326
  )

  expect_error(
    predict_axis_surface(model_fit, grid_lonlat),
    "longitude/latitude"
  )
})
