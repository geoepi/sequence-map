test_that("alignment_to_variant_matrix drops invariant sites by default", {
  alignment <- ape::as.DNAbin(matrix(
    c(
      "a", "c", "g",
      "a", "t", "g",
      "a", "c", "g"
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("seq1", "seq2", "seq3"), NULL)
  ))

  result <- alignment_to_variant_matrix(alignment)

  expect_equal(nrow(result$variant_matrix), 3)
  expect_equal(ncol(result$variant_matrix), 4)
  expect_equal(
    colnames(result$variant_matrix),
    c("site_2_A", "site_2_C", "site_2_G", "site_2_T")
  )
  expect_true(result$site_summary$retained[2])
  expect_false(result$site_summary$retained[1])
  expect_false(result$site_summary$retained[3])
  expect_equal(result$site_summary$site_class, c("invariant_canonical", "variable_canonical", "invariant_canonical"))
})

test_that("alignment_to_variant_matrix tracks ambiguous and gap states", {
  alignment <- ape::as.DNAbin(matrix(
    c(
      "a", "-", "n",
      "c", "r", "t"
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("seq1", "seq2"), NULL)
  ))

  result <- alignment_to_variant_matrix(alignment, drop_invariant = FALSE)

  expect_equal(result$site_summary$n_gap, c(0, 1, 0))
  expect_equal(result$site_summary$n_ambiguous, c(0, 1, 0))
  expect_equal(result$site_summary$n_missing, c(0, 0, 1))
  expect_equal(result$sequence_summary$n_gap, c(1, 0))
  expect_equal(result$sequence_summary$n_ambiguous, c(0, 1))
  expect_equal(result$sequence_summary$n_missing, c(1, 0))
  expect_false(any(grepl("^site_2_", colnames(result$variant_matrix))))
  expect_true(all(grepl("^site_[13]_", colnames(result$variant_matrix))))
  expect_equal(result$site_summary$site_class, c("variable_canonical", "noncanonical_only", "invariant_canonical_with_noncanonical"))
})

test_that("alignment_to_variant_matrix requires DNAbin input", {
  expect_error(
    alignment_to_variant_matrix(matrix(c("A", "C"), nrow = 1)),
    "DNAbin"
  )
})

test_that("alignment_to_variant_matrix retains expected variable sites", {
  alignment <- ape::as.DNAbin(matrix(
    c(
      "a", "c", "-", "t",
      "a", "g", "n", "t",
      "a", "c", "-", "a"
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("seq1", "seq2", "seq3"), NULL)
  ))

  result <- alignment_to_variant_matrix(alignment, drop_invariant = TRUE)

  expect_equal(nrow(result$variant_matrix), 3)
  expect_equal(ncol(result$variant_matrix), 8)
  expect_equal(
    colnames(result$variant_matrix),
    c(
      "site_2_A", "site_2_C", "site_2_G", "site_2_T",
      "site_4_A", "site_4_C", "site_4_G", "site_4_T"
    )
  )
  expect_equal(which(result$site_summary$retained), c(2, 4))
  expect_equal(
    result$site_summary$site_class,
    c(
      "invariant_canonical",
      "variable_canonical",
      "noncanonical_only",
      "variable_canonical"
    )
  )
})

test_that("alignment_to_variant_matrix keeps invariant canonical sites when requested", {
  alignment <- ape::as.DNAbin(matrix(
    c(
      "a", "c", "-", "t",
      "a", "g", "n", "t",
      "a", "c", "-", "a"
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("seq1", "seq2", "seq3"), NULL)
  ))

  result <- alignment_to_variant_matrix(alignment, drop_invariant = FALSE)

  expect_equal(which(result$site_summary$retained), c(1, 2, 4))
  expect_equal(ncol(result$variant_matrix), 12)
  expect_false(result$site_summary$retained[3])
})

test_that("alignment_to_variant_matrix errors when no sites remain", {
  alignment <- ape::as.DNAbin(matrix(
    c(
      "-", "n", "a",
      "-", "n", "a",
      "-", "n", "a"
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("seq1", "seq2", "seq3"), NULL)
  ))

  expect_error(
    alignment_to_variant_matrix(alignment),
    "No retained sites remain after filtering"
  )
})
