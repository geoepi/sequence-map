test_that("alignment IDs must be present in metadata", {
  alignment <- list(a = c("A", "C"), b = c("G", "T"), c = c("T", "T"))

  metadata <- data.frame(
    sequence_id = c("a", "b"),
    longitude = c(1, 2),
    latitude = c(3, 4),
    stringsAsFactors = FALSE
  )

  expect_error(
    validate_alignment_metadata_match(alignment, metadata),
    "c"
  )
})

test_that("coordinate columns must be numeric", {
  metadata <- data.frame(
    sequence_id = c("a", "b"),
    longitude = c("10.5", "west"),
    latitude = c("5.2", "6.3"),
    stringsAsFactors = FALSE
  )

  expect_error(
    validate_sequence_metadata(metadata),
    "Coordinate columns must be numeric"
  )
})

test_that("retained sequences must have complete coordinates", {
  alignment <- list(a = c("A", "C"), b = c("G", "T"))

  metadata <- data.frame(
    sequence_id = c("a", "b", "c"),
    longitude = c(10, NA, 30),
    latitude = c(5, 6, NA),
    stringsAsFactors = FALSE
  )

  expect_error(
    validate_sequence_metadata(metadata, alignment = alignment),
    "b"
  )
})
