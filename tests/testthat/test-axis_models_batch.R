test_that("batch axis model helpers skip gracefully if INLA is not installed", {
  if (requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is installed; skip not-installed branch.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(PC1_mean = c(0, 1, 2), x = c(0, 100, 200), y = c(0, 100, 0)),
    coords = c("x", "y"),
    crs = 3857
  )

  expect_error(
    fit_axis_models(points_sf = points_sf, mesh = list(mesh = list(loc = matrix(0, nrow = 3, ncol = 2)))),
    "Package 'INLA' is required"
  )
})

test_that("fit_axis_models auto-detects PC and LD mean columns and excludes summary columns", {
  data <- sf::st_as_sf(
    data.frame(
      PC1_mean = c(0.1, 0.3, 0.8, 1),
      PC2_mean = c(0.2, 0.4, 0.7, 0.9),
      LD1_mean = c(1.1, 1.3, 1.8, 2),
      PC1_sd = c(0.1, 0.1, 0.1, 0.1),
      LD1_pred = c(1, 1, 1, 1),
      mean = c(0, 0, 0, 0),
      x = c(0, 1000, 0, 1000),
      y = c(0, 0, 1000, 1000)
    ),
    coords = c("x", "y"),
    crs = 3857
  )

  detected <- detect_axis_mean_response_cols(data)

  expect_equal(detected, c("PC1_mean", "PC2_mean", "LD1_mean"))
})

test_that("fit_axis_models returns model objects and model summaries", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(
      PC1_mean = c(0.1, 0.3, 0.8, 1.0),
      LD1_mean = c(1.1, 1.3, 1.8, 2.0),
      x = c(0, 1000, 0, 1000),
      y = c(0, 0, 1000, 1000)
    ),
    coords = c("x", "y"),
    crs = 3857
  )
  mesh_obj <- make_inla_mesh(points_sf, max_edge = c(800, 1600), cutoff = 1, offset = c(200, 400), crs_projected = 3857)

  result <- fit_axis_models(points_sf = points_sf, mesh = mesh_obj)

  expect_s3_class(result, "sequencemap_axis_models")
  expect_equal(sort(result$response_cols), c("LD1_mean", "PC1_mean"))
  expect_equal(result$n_successful, 2)
  expect_true(all(c("response_col", "n_observations", "n_mesh_vertices", "waic", "dic") %in% names(result$model_summaries)))
  expect_true(all(result$successful_response_cols %in% names(result$models)))
})

test_that("fit_axis_models records failures when continue_on_error is true", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(
      PC1_mean = c(0.1, 0.3, 0.8, 1.0),
      PC2_mean = c(NA, NA, 0.2, NA),
      x = c(0, 1000, 0, 1000),
      y = c(0, 0, 1000, 1000)
    ),
    coords = c("x", "y"),
    crs = 3857
  )
  mesh_obj <- make_inla_mesh(
    points_sf = sf::st_as_sf(
      data.frame(PC1_mean = c(0.1, 0.2, 0.3, 0.4), x = c(0, 1000, 0, 1000), y = c(0, 0, 1000, 1000)),
      coords = c("x", "y"),
      crs = 3857
    ),
    max_edge = c(800, 1600),
    cutoff = 1,
    offset = c(200, 400),
    crs_projected = 3857
  )

  expect_warning(
    result <- fit_axis_models(
      points_sf = points_sf,
      response_cols = c("PC1_mean", "PC2_mean"),
      mesh = mesh_obj,
      continue_on_error = TRUE
    ),
    "Dropped 3 observations"
  )

  expect_equal(result$n_successful, 1)
  expect_equal(result$n_failed, 1)
  expect_true("PC2_mean" %in% names(result$errors))
})

test_that("fit_axis_models stops on failures when continue_on_error is false", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(
      PC1_mean = c(0.1, 0.3, 0.8, 1.0),
      PC2_mean = c(NA, NA, 0.2, NA),
      x = c(0, 1000, 0, 1000),
      y = c(0, 0, 1000, 1000)
    ),
    coords = c("x", "y"),
    crs = 3857
  )
  mesh_obj <- make_inla_mesh(
    points_sf = sf::st_as_sf(
      data.frame(PC1_mean = c(0.1, 0.2, 0.3, 0.4), x = c(0, 1000, 0, 1000), y = c(0, 0, 1000, 1000)),
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
      fit_axis_models(
        points_sf = points_sf,
        response_cols = c("PC1_mean", "PC2_mean"),
        mesh = mesh_obj,
        continue_on_error = FALSE
      ),
      "At least 3 non-missing observations"
    ),
    "Dropped 3 observations"
  )
})

test_that("predict_axis_surfaces returns one surface per model and a combined surface", {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    skip("INLA is not installed.")
  }

  points_sf <- sf::st_as_sf(
    data.frame(
      PC1_mean = c(0.1, 0.3, 0.8, 1.0),
      LD1_mean = c(1.1, 1.3, 1.8, 2.0),
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
  mesh_obj <- make_inla_mesh(points_sf, boundary = boundary, prediction_grid = grid, max_edge = c(800, 1600), cutoff = 1, offset = c(200, 400), crs_projected = 3857)
  model_batch <- fit_axis_models(points_sf = points_sf, mesh = mesh_obj)

  surface_batch <- predict_axis_surfaces(model_batch, grid)

  expect_s3_class(surface_batch, "sequencemap_axis_surfaces")
  expect_equal(length(surface_batch$surfaces), 2)
  expect_true(inherits(surface_batch$combined_surface, "sf"))
  expect_true(all(c("grid_id", "x", "y", "geometry", "PC1_mean_pred", "PC1_mean_sd", "LD1_mean_pred", "LD1_mean_sd") %in% names(surface_batch$combined_surface)))
  expect_false(any(c("mean", "sd", "q025", "q975") %in% names(surface_batch$combined_surface)))
})

test_that("predict_axis_surfaces records prediction failures when continue_on_error is true", {
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
  mesh_obj <- make_inla_mesh(points_sf, boundary = boundary, prediction_grid = grid, max_edge = c(800, 1600), cutoff = 1, offset = c(200, 400), crs_projected = 3857)
  valid_model <- fit_axis_model(points_sf, "PC1_mean", mesh = mesh_obj)
  invalid_model <- list(response_col = "LD1_mean")

  surface_batch <- predict_axis_surfaces(
    axis_models = list(PC1_mean = valid_model, LD1_mean = invalid_model),
    prediction_grid = grid,
    continue_on_error = TRUE
  )

  expect_equal(surface_batch$successful_response_cols, "PC1_mean")
  expect_equal(surface_batch$failed_response_cols, "LD1_mean")
  expect_true("LD1_mean" %in% names(surface_batch$errors))
})
