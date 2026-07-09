make_polygon_axis_summaries <- function() {
  data.frame(
    polygon_row = c(1L, 2L),
    polygon_id = c("province_a", "province_b"),
    n_grid_cells = c(20L, 15L),
    PC1_mean_pred_mean = c(3, 6),
    PC2_mean_pred_mean = c(4, NA_real_),
    PC3_mean_pred_mean = c(12, 8),
    PC1_mean_pred_sd = c(1, 2),
    PC2_mean_pred_sd = c(2, NA_real_),
    PC3_mean_pred_sd = c(2, 3),
    PC1_mean_sd_mean = c(0.1, 0.4),
    PC2_mean_sd_mean = c(0.2, NA_real_),
    PC3_mean_sd_mean = c(0.3, 0.6),
    stringsAsFactors = FALSE
  )
}

test_that("calculate_polygon_axis_metrics detects PC summary axes", {
  result <- calculate_polygon_axis_metrics(make_polygon_axis_summaries())

  expect_equal(result$n_axes_used, c(3, 2))
  expect_true(all(c(
    "axis_centroid_distance",
    "axis_surface_dispersion",
    "mean_axis_uncertainty",
    "n_axes_used"
  ) %in% names(result)))
})

test_that("calculate_polygon_axis_metrics calculates centroid distance", {
  result <- calculate_polygon_axis_metrics(make_polygon_axis_summaries())

  expect_equal(result$axis_centroid_distance, c(13, 10))
})

test_that("calculate_polygon_axis_metrics calculates surface dispersion", {
  result <- calculate_polygon_axis_metrics(make_polygon_axis_summaries())

  expect_equal(result$axis_surface_dispersion, c(9, 13))
})

test_that("calculate_polygon_axis_metrics calculates mean uncertainty", {
  result <- calculate_polygon_axis_metrics(make_polygon_axis_summaries())

  expect_equal(result$mean_axis_uncertainty, c(0.2, 0.5))
})

test_that("calculate_polygon_axis_metrics handles missing values gracefully", {
  summaries <- make_polygon_axis_summaries()
  summaries[2, c(
    "PC1_mean_pred_mean",
    "PC3_mean_pred_mean",
    "PC1_mean_pred_sd",
    "PC3_mean_pred_sd",
    "PC1_mean_sd_mean",
    "PC3_mean_sd_mean"
  )] <- NA_real_

  result <- calculate_polygon_axis_metrics(summaries)

  expect_equal(result$n_axes_used[2], 0)
  expect_true(is.na(result$axis_centroid_distance[2]))
  expect_true(is.na(result$axis_surface_dispersion[2]))
  expect_true(is.na(result$mean_axis_uncertainty[2]))
})

test_that("calculate_polygon_axis_metrics errors without usable axes", {
  expect_error(
    calculate_polygon_axis_metrics(data.frame(polygon_id = "province_a")),
    "No usable predicted axis mean columns"
  )
})

test_that("calculate_polygon_axis_metrics preserves polygon identifiers and cells", {
  summaries <- make_polygon_axis_summaries()
  result <- calculate_polygon_axis_metrics(summaries)

  expect_equal(result$polygon_id, summaries$polygon_id)
  expect_equal(result$n_grid_cells, summaries$n_grid_cells)
  expect_equal(result$PC1_mean_pred_mean, summaries$PC1_mean_pred_mean)
})
