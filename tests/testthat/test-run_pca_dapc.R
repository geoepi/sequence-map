test_that("run_sequence_pca returns tidy scores and variance summaries", {
  encoded <- list(
    variant_matrix = matrix(
      c(
        1, 0, 0, 1,
        1, 0, 1, 0,
        0, 1, 1, 0
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("seq1", "seq2", "seq3"), paste0("v", 1:4))
    )
  )
  metadata <- data.frame(
    sequence_id = c("seq1", "seq2", "seq3"),
    region = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )

  result <- run_sequence_pca(encoded, metadata = metadata)

  expect_true(all(c("sequence_id", "PC1") %in% names(result$scores)))
  expect_equal(result$scores$sequence_id, c("seq1", "seq2", "seq3"))
  expect_true("region" %in% names(result$scores))
  expect_equal(sum(result$variance_summary$proportion_variance), 1, tolerance = 1e-8)
})

test_that("run_sequence_dapc uses supplied grouping metadata", {
  encoded <- list(
    variant_matrix = matrix(
      c(
        1, 0, 1, 0,
        1, 0, 0.9, 0.1,
        0, 1, 0, 1,
        0.1, 0.9, 0, 1
      ),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(c("seq1", "seq2", "seq3", "seq4"), paste0("v", 1:4))
    )
  )
  metadata <- data.frame(
    sequence_id = c("seq1", "seq2", "seq3", "seq4"),
    cluster = c("alpha", "alpha", "beta", "beta"),
    stringsAsFactors = FALSE
  )

  result <- run_sequence_dapc(
    encoded,
    metadata = metadata,
    grouping_col = "cluster",
    n_pca = 2,
    n_da = 1
  )

  expect_true(all(c("sequence_id", "group", "assigned_group", "LD1") %in% names(result$scores)))
  expect_equal(sort(unique(result$scores$group)), c("alpha", "beta"))
  expect_true(all(result$eigen_summary$eigenvalue >= 0))
  expect_false(result$inferred_groups)
})

test_that("run_sequence_dapc can infer groups", {
  encoded <- list(
    variant_matrix = matrix(
      c(
        1, 0, 1, 0,
        0.9, 0.1, 1, 0,
        0, 1, 0, 1,
        0, 1, 0.1, 0.9
      ),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(c("seq1", "seq2", "seq3", "seq4"), paste0("v", 1:4))
    )
  )
  result <- run_sequence_dapc(
    encoded,
    infer_groups = TRUE,
    n_pca = 2,
    n_da = 1,
    n_clust = 2,
    seed = 1
  )

  expect_true(result$inferred_groups)
  expect_equal(length(unique(result$scores$group)), 2)
  expect_true("LD1" %in% names(result$scores))
})

test_that("run_sequence_pca errors for zero-variance matrices", {
  encoded <- list(
    variant_matrix = matrix(
      c(
        1, 0, 1, 0,
        1, 0, 1, 0,
        1, 0, 1, 0
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("seq1", "seq2", "seq3"), paste0("v", 1:4))
    )
  )

  expect_error(
    run_sequence_pca(encoded),
    "nonzero variance"
  )
})

test_that("run_sequence_dapc errors for empty feature matrices", {
  encoded <- list(
    variant_matrix = matrix(
      numeric(0),
      nrow = 2,
      ncol = 0,
      dimnames = list(c("seq1", "seq2"), NULL)
    )
  )

  expect_error(
    run_sequence_dapc(encoded, infer_groups = TRUE),
    "at least one encoded feature"
  )
})
