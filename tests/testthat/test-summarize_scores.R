test_that("summarize_location_scores averages PC and LD scores by location", {
  scores <- data.frame(
    sequence_id = c("s1", "s2", "s3"),
    PC1 = c(1, 3, 5),
    LD1 = c(2, 4, 6),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    sequence_id = c("s1", "s2", "s3"),
    location = c("loc1", "loc1", "loc2"),
    longitude = c(10, 10, 20),
    latitude = c(50, 50, 60),
    stringsAsFactors = FALSE
  )

  result <- summarize_location_scores(
    scores,
    metadata = metadata,
    location_col = "location"
  )

  expect_equal(result$n_sequences, c(2, 1))
  expect_equal(result$PC1_mean, c(2, 5))
  expect_equal(result$LD1_mean, c(3, 6))
})

test_that("summarize_location_scores groups by coordinates when no location column is supplied", {
  scores <- data.frame(
    sequence_id = c("s1", "s2", "s3"),
    PC1 = c(1, 3, 5),
    longitude = c(10, 10, 20),
    latitude = c(50, 50, 60),
    stringsAsFactors = FALSE
  )

  result <- summarize_location_scores(scores)

  expect_equal(nrow(result), 2)
  expect_equal(result$n_sequences, c(2, 1))
  expect_equal(result$PC1_mean, c(2, 5))
})

test_that("summarize_location_scores returns NA sd for single-sequence locations", {
  scores <- data.frame(
    sequence_id = c("s1", "s2"),
    PC1 = c(1, 2),
    longitude = c(10, 20),
    latitude = c(50, 60),
    stringsAsFactors = FALSE
  )

  result <- summarize_location_scores(scores)

  expect_true(all(is.na(result$PC1_sd)))
})

test_that("summarize_location_diversity returns div_trace for locations with enough sequences", {
  scores <- data.frame(
    sequence_id = c("s1", "s2", "s3"),
    PC1 = c(0, 2, 5),
    PC2 = c(0, 2, 5),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    sequence_id = c("s1", "s2", "s3"),
    location = c("loc1", "loc1", "loc2"),
    longitude = c(10, 10, 20),
    latitude = c(50, 50, 60),
    stringsAsFactors = FALSE
  )

  result <- summarize_location_diversity(
    scores,
    metadata = metadata,
    location_col = "location",
    transform = "none"
  )

  expect_equal(result$n_sequences, c(2, 1))
  expect_true(!is.na(result$div_trace[1]))
  expect_true(result$div_trace[1] > 0)
})

test_that("summarize_location_diversity returns NA when too few sequences are present", {
  scores <- data.frame(
    sequence_id = c("s1", "s2"),
    PC1 = c(0, 5),
    PC2 = c(0, 5),
    longitude = c(10, 20),
    latitude = c(50, 60),
    stringsAsFactors = FALSE
  )

  result <- summarize_location_diversity(scores, min_sequences = 2)

  expect_true(all(is.na(result$div_trace)))
  expect_true(all(is.na(result$mean_pairwise_ordination_distance)))
})

test_that("summarize_location_diversity auto-detects PC columns only", {
  scores <- data.frame(
    sequence_id = c("s1", "s2"),
    PC1 = c(0, 1),
    PC_2 = c(1, 0),
    LD1 = c(5, 6),
    posterior_alpha = c(0.4, 0.6),
    region = c("x", "y"),
    longitude = c(10, 10),
    latitude = c(50, 50),
    stringsAsFactors = FALSE
  )

  result <- summarize_location_diversity(scores, min_sequences = 2)

  expect_equal(result$n_axes, 2)
})

test_that("summary functions error clearly for missing IDs, coordinates, and score columns", {
  scores <- data.frame(
    sequence_id = c("s1", "s2"),
    PC1 = c(1, 2),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    sequence_id = c("s1"),
    location = "loc1",
    longitude = 10,
    latitude = 50,
    stringsAsFactors = FALSE
  )

  expect_error(
    summarize_location_scores(scores, metadata = metadata, location_col = "location"),
    "missing from metadata"
  )

  expect_error(
    summarize_location_scores(scores),
    "Coordinates are required"
  )

  expect_error(
    summarize_location_diversity(
      data.frame(
        sequence_id = c("s1", "s2"),
        posterior_alpha = c(0.4, 0.6),
        longitude = c(10, 10),
        latitude = c(50, 50),
        stringsAsFactors = FALSE
      )
    ),
    "No valid score columns"
  )
})
