#' Validate sequence metadata columns
#'
#' Confirms that required identifier and coordinate columns exist, that
#' coordinates are numeric, and that retained sequences have complete
#' coordinates.
#'
#' @param metadata A data frame containing sequence metadata.
#' @param alignment Optional alignment object returned by [read_alignment()].
#' @param id_col Name of the metadata column containing sequence identifiers.
#' @param longitude_col Name of the metadata longitude column.
#' @param latitude_col Name of the metadata latitude column.
#'
#' @return The input metadata, invisibly.
#' @export
validate_sequence_metadata <- function(
  metadata,
  alignment = NULL,
  id_col = "sequence_id",
  longitude_col = "longitude",
  latitude_col = "latitude"
) {
  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data frame.", call. = FALSE)
  }

  required_columns <- c(id_col, longitude_col, latitude_col)
  missing_columns <- setdiff(required_columns, names(metadata))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "Metadata is missing required columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  longitude_numeric <- suppressWarnings(as.numeric(metadata[[longitude_col]]))
  latitude_numeric <- suppressWarnings(as.numeric(metadata[[latitude_col]]))

  invalid_longitude <- !is.na(metadata[[longitude_col]]) & is.na(longitude_numeric)
  invalid_latitude <- !is.na(metadata[[latitude_col]]) & is.na(latitude_numeric)

  if (any(invalid_longitude) || any(invalid_latitude)) {
    invalid_ids <- unique(c(
      metadata[[id_col]][invalid_longitude],
      metadata[[id_col]][invalid_latitude]
    ))
    stop(
      sprintf(
        "Coordinate columns must be numeric. Offending sequence IDs: %s",
        format_example_ids(invalid_ids)
      ),
      call. = FALSE
    )
  }

  retained_metadata <- metadata
  if (!is.null(alignment)) {
    validate_alignment_metadata_match(alignment = alignment, metadata = metadata, id_col = id_col)
    retained_ids <- names(alignment)
    retained_metadata <- metadata[metadata[[id_col]] %in% retained_ids, , drop = FALSE]
  }

  retained_longitude <- suppressWarnings(as.numeric(retained_metadata[[longitude_col]]))
  retained_latitude <- suppressWarnings(as.numeric(retained_metadata[[latitude_col]]))
  missing_coordinates <- is.na(retained_longitude) | is.na(retained_latitude)
  retained_missing <- retained_metadata[[id_col]][missing_coordinates]

  if (length(retained_missing) > 0L) {
    stop(
      sprintf(
        "Retained sequences must have complete coordinates. Missing coordinates for: %s",
        format_example_ids(unique(retained_missing))
      ),
      call. = FALSE
    )
  }

  invisible(metadata)
}

#' Validate alignment and metadata identifier matching
#'
#' Confirms that every sequence ID in an alignment is present in the metadata.
#'
#' @param alignment Alignment object returned by [read_alignment()].
#' @param metadata A data frame containing sequence metadata.
#' @param id_col Name of the metadata column containing sequence identifiers.
#'
#' @return `TRUE`, invisibly.
#' @export
validate_alignment_metadata_match <- function(
  alignment,
  metadata,
  id_col = "sequence_id"
) {
  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data frame.", call. = FALSE)
  }

  if (!id_col %in% names(metadata)) {
    stop(
      sprintf("Metadata is missing the sequence ID column '%s'.", id_col),
      call. = FALSE
    )
  }

  sequence_ids <- names(alignment)
  if (is.null(sequence_ids) || length(sequence_ids) == 0L) {
    stop("`alignment` must contain named sequences.", call. = FALSE)
  }

  metadata_ids <- metadata[[id_col]]
  missing_ids <- setdiff(sequence_ids, metadata_ids)

  if (length(missing_ids) > 0L) {
    stop(
      sprintf(
        "Sequence IDs in the alignment are missing from metadata (%d missing). Examples: %s",
        length(missing_ids),
        format_example_ids(missing_ids)
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Format ID examples for validation messages
#'
#' @param ids Character vector of sequence identifiers.
#'
#' @return A length-1 character string.
#' @keywords internal
format_example_ids <- function(ids) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]

  if (length(ids) == 0L) {
    return("<none>")
  }

  shown_ids <- utils::head(ids, 20L)
  suffix <- if (length(ids) > 20L) ", ..." else ""
  paste0(paste(shown_ids, collapse = ", "), suffix)
}
