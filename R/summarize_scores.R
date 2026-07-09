#' Summarize ordination scores by sampled location
#'
#' Aggregates sequence-level PCA and/or DAPC scores to unique sampled
#' locations. These are summary statistics of ordination scores, not spatial
#' model outputs.
#'
#' @param scores A data frame containing `sequence_id` and one or more score
#'   columns.
#' @param metadata Optional metadata data frame containing `sequence_id` and
#'   location or coordinate columns.
#' @param sequence_id_col Name of the sequence identifier column.
#' @param location_col Optional location identifier column.
#' @param lon_col Name of the longitude column.
#' @param lat_col Name of the latitude column.
#' @param score_cols Optional score columns to summarize. If `NULL`, columns
#'   beginning with `PC`, `PC_`, `LD`, or `LD_` are detected automatically.
#' @param extra_summary_cols Optional metadata columns to carry forward as
#'   location-level summaries when present.
#'
#' @return A tidy data frame with one row per sampled location.
#' @export
summarize_location_scores <- function(
  scores,
  metadata = NULL,
  sequence_id_col = "sequence_id",
  location_col = NULL,
  lon_col = "longitude",
  lat_col = "latitude",
  score_cols = NULL,
  extra_summary_cols = NULL
) {
  prepared <- prepare_score_summary_data(
    scores = scores,
    metadata = metadata,
    sequence_id_col = sequence_id_col,
    location_col = location_col,
    lon_col = lon_col,
    lat_col = lat_col,
    score_cols = score_cols,
    score_mode = "all"
  )

  score_data <- prepared$data
  split_groups <- split(score_data, score_data$.location_group, drop = TRUE)

  rows <- lapply(split_groups, function(group_data) {
    location_info <- summarize_location_coordinates(
      group_data = group_data,
      location_col = location_col,
      lon_col = lon_col,
      lat_col = lat_col
    )

    score_summary <- lapply(prepared$score_cols, function(score_col) {
      values <- group_data[[score_col]]
      c(
        mean = mean(values, na.rm = TRUE),
        sd = if (length(values) >= 2L) stats::sd(values, na.rm = TRUE) else NA_real_,
        min = min(values, na.rm = TRUE),
        max = max(values, na.rm = TRUE)
      )
    })
    names(score_summary) <- prepared$score_cols

    row <- c(
      location_info,
      list(n_sequences = nrow(group_data))
    )

    for (score_col in prepared$score_cols) {
      row[[paste0(score_col, "_mean")]] <- unname(score_summary[[score_col]]["mean"])
      row[[paste0(score_col, "_sd")]] <- unname(score_summary[[score_col]]["sd"])
      row[[paste0(score_col, "_min")]] <- unname(score_summary[[score_col]]["min"])
      row[[paste0(score_col, "_max")]] <- unname(score_summary[[score_col]]["max"])
    }

    if (!is.null(extra_summary_cols)) {
      available_extra <- extra_summary_cols[extra_summary_cols %in% names(group_data)]
      for (column_name in available_extra) {
        row[[column_name]] <- summarize_extra_column(group_data[[column_name]])
      }
    }

    as.data.frame(row, stringsAsFactors = FALSE)
  })

  do.call(rbind, rows)
}

#' Summarize ordination-space diversity by sampled location
#'
#' Calculates within-location dispersion in reduced ordination space. This is
#' an ordination-space diversity summary, not nucleotide diversity.
#'
#' @param scores A sequence-level ordination score table.
#' @param metadata Optional metadata data frame containing `sequence_id` and
#'   location or coordinate columns.
#' @param sequence_id_col Name of the sequence identifier column.
#' @param location_col Optional location identifier column.
#' @param lon_col Name of the longitude column.
#' @param lat_col Name of the latitude column.
#' @param score_cols Optional score columns. If `NULL`, principal component
#'   columns are detected automatically.
#' @param min_sequences Minimum sequences required to compute dispersion.
#' @param transform Transformation applied to `div_trace`.
#'
#' @return A tidy data frame with one row per sampled location.
#' @export
summarize_location_diversity <- function(
  scores,
  metadata = NULL,
  sequence_id_col = "sequence_id",
  location_col = NULL,
  lon_col = "longitude",
  lat_col = "latitude",
  score_cols = NULL,
  min_sequences = 2,
  transform = c("none", "log", "log1p")
) {
  transform <- match.arg(transform)

  prepared <- prepare_score_summary_data(
    scores = scores,
    metadata = metadata,
    sequence_id_col = sequence_id_col,
    location_col = location_col,
    lon_col = lon_col,
    lat_col = lat_col,
    score_cols = score_cols,
    score_mode = "pca_only"
  )

  score_data <- prepared$data
  split_groups <- split(score_data, score_data$.location_group, drop = TRUE)

  rows <- lapply(split_groups, function(group_data) {
    location_info <- summarize_location_coordinates(
      group_data = group_data,
      location_col = location_col,
      lon_col = lon_col,
      lat_col = lat_col
    )

    n_sequences <- nrow(group_data)
    n_axes <- length(prepared$score_cols)

    row <- c(
      location_info,
      list(
        n_sequences = n_sequences,
        n_axes = n_axes
      )
    )

    if (n_sequences < min_sequences) {
      row$div_trace <- NA_real_
      row$mean_pairwise_ordination_distance <- NA_real_
      row$mean_distance_to_location_centroid <- NA_real_
      row$transformed_div_trace <- NA_real_
      return(as.data.frame(row, stringsAsFactors = FALSE))
    }

    score_matrix <- as.matrix(group_data[, prepared$score_cols, drop = FALSE])
    storage.mode(score_matrix) <- "double"

    covariance_matrix <- stats::cov(score_matrix)
    div_trace <- sum(diag(covariance_matrix))

    if (n_sequences >= 2L) {
      pairwise_distances <- stats::dist(score_matrix)
      mean_pairwise_distance <- if (length(pairwise_distances) > 0L) {
        mean(as.numeric(pairwise_distances))
      } else {
        0
      }
    } else {
      mean_pairwise_distance <- NA_real_
    }

    centroid <- colMeans(score_matrix)
    centered <- sweep(score_matrix, 2, centroid, FUN = "-")
    mean_distance_to_centroid <- mean(sqrt(rowSums(centered^2)))

    row$div_trace <- div_trace
    row$mean_pairwise_ordination_distance <- mean_pairwise_distance
    row$mean_distance_to_location_centroid <- mean_distance_to_centroid
    row$transformed_div_trace <- transform_div_trace(div_trace, transform)

    as.data.frame(row, stringsAsFactors = FALSE)
  })

  do.call(rbind, rows)
}

#' Prepare ordination score data for location summaries
#'
#' @param scores Score data frame.
#' @param metadata Optional metadata data frame.
#' @param sequence_id_col Sequence identifier column name.
#' @param location_col Optional location column.
#' @param lon_col Longitude column.
#' @param lat_col Latitude column.
#' @param score_cols Optional explicit score columns.
#' @param score_mode Score autodetection mode.
#'
#' @return A list with joined `data` and selected `score_cols`.
#' @keywords internal
prepare_score_summary_data <- function(
  scores,
  metadata,
  sequence_id_col,
  location_col,
  lon_col,
  lat_col,
  score_cols,
  score_mode = c("all", "pca_only")
) {
  score_mode <- match.arg(score_mode)

  if (!is.data.frame(scores)) {
    stop("`scores` must be a data frame.", call. = FALSE)
  }

  if (!sequence_id_col %in% names(scores)) {
    stop(
      sprintf("`scores` is missing the sequence ID column '%s'.", sequence_id_col),
      call. = FALSE
    )
  }

  data <- scores
  if (sequence_id_col != "sequence_id") {
    names(data)[names(data) == sequence_id_col] <- "sequence_id"
  }

  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) {
      stop("`metadata` must be a data frame when supplied.", call. = FALSE)
    }

    if (!sequence_id_col %in% names(metadata)) {
      stop(
        sprintf("`metadata` is missing the sequence ID column '%s'.", sequence_id_col),
        call. = FALSE
      )
    }

    metadata_copy <- metadata
    names(metadata_copy)[names(metadata_copy) == sequence_id_col] <- "sequence_id"

    match_index <- match(data$sequence_id, metadata_copy$sequence_id)
    missing_ids <- data$sequence_id[is.na(match_index)]
    if (length(missing_ids) > 0L) {
      stop(
        sprintf(
          "Sequence IDs in `scores` are missing from metadata (%d missing). Examples: %s",
          length(missing_ids),
          format_example_ids(missing_ids)
        ),
        call. = FALSE
      )
    }

    metadata_aligned <- metadata_copy[match_index, , drop = FALSE]
    metadata_aligned <- metadata_aligned[, setdiff(names(metadata_aligned), "sequence_id"), drop = FALSE]
    duplicated_columns <- intersect(names(data), names(metadata_aligned))
    if (length(duplicated_columns) > 0L) {
      metadata_aligned <- metadata_aligned[, setdiff(names(metadata_aligned), duplicated_columns), drop = FALSE]
    }
    data <- cbind(data, metadata_aligned)
  }

  selected_score_cols <- detect_score_columns(
    data = data,
    score_cols = score_cols,
    score_mode = score_mode
  )

  validate_coordinate_columns(
    data = data,
    location_col = location_col,
    lon_col = lon_col,
    lat_col = lat_col
  )

  data[[lon_col]] <- coerce_numeric_column(data[[lon_col]], lon_col)
  data[[lat_col]] <- coerce_numeric_column(data[[lat_col]], lat_col)

  if (!is.null(location_col)) {
    if (!location_col %in% names(data)) {
      stop(
        sprintf("Location column '%s' was not found after joining metadata.", location_col),
        call. = FALSE
      )
    }
    data$.location_group <- as.character(data[[location_col]])
  } else {
    data$.location_group <- paste(data[[lon_col]], data[[lat_col]], sep = "||")
  }

  data$.location_group[is.na(data$.location_group) | !nzchar(data$.location_group)] <- NA_character_
  if (anyNA(data$.location_group)) {
    bad_ids <- data$sequence_id[is.na(data$.location_group)]
    stop(
      sprintf(
        "Location grouping is missing for some sequences. Examples: %s",
        format_example_ids(bad_ids)
      ),
      call. = FALSE
    )
  }

  list(data = data, score_cols = selected_score_cols)
}

#' Detect score columns from a score table
#'
#' @param data A data frame.
#' @param score_cols Optional explicit columns.
#' @param score_mode Autodetection mode.
#'
#' @return Character vector of score columns.
#' @keywords internal
detect_score_columns <- function(data, score_cols = NULL, score_mode = c("all", "pca_only")) {
  score_mode <- match.arg(score_mode)

  if (!is.null(score_cols)) {
    missing_cols <- setdiff(score_cols, names(data))
    if (length(missing_cols) > 0L) {
      stop(
        sprintf("Requested score columns were not found: %s", paste(missing_cols, collapse = ", ")),
        call. = FALSE
      )
    }
    selected <- score_cols
  } else {
    patterns <- if (score_mode == "all") {
      c("^PC_?[0-9]+$", "^LD_?[0-9]+$")
    } else {
      c("^PC_?[0-9]+$")
    }

    selected <- unique(unlist(lapply(patterns, function(pattern) {
      grep(pattern, names(data), value = TRUE)
    })))
  }

  selected <- selected[vapply(data[selected], is.numeric, logical(1))]
  if (length(selected) == 0L) {
    stop("No valid score columns were found for summarization.", call. = FALSE)
  }

  selected
}

#' Validate required coordinate columns
#'
#' @param data Joined score data.
#' @param location_col Optional location column.
#' @param lon_col Longitude column name.
#' @param lat_col Latitude column name.
#'
#' @return `TRUE`, invisibly.
#' @keywords internal
validate_coordinate_columns <- function(data, location_col, lon_col, lat_col) {
  missing_cols <- setdiff(c(lon_col, lat_col), names(data))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "Coordinates are required but missing columns were not found: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(location_col) && !location_col %in% names(data)) {
    stop(
      sprintf("Location column '%s' was not found after joining metadata.", location_col),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Coerce a coordinate column to numeric
#'
#' @param x Column vector.
#' @param column_name Column name.
#'
#' @return Numeric vector.
#' @keywords internal
coerce_numeric_column <- function(x, column_name) {
  numeric_x <- suppressWarnings(as.numeric(x))
  invalid <- !is.na(x) & is.na(numeric_x)
  if (any(invalid)) {
    stop(
      sprintf("Coordinate column '%s' must be numeric.", column_name),
      call. = FALSE
    )
  }
  numeric_x
}

#' Summarize location coordinates and identifiers
#'
#' @param group_data Grouped data frame.
#' @param location_col Optional location column.
#' @param lon_col Longitude column.
#' @param lat_col Latitude column.
#'
#' @return Named list.
#' @keywords internal
summarize_location_coordinates <- function(group_data, location_col, lon_col, lat_col) {
  longitude <- group_data[[lon_col]]
  latitude <- group_data[[lat_col]]

  row <- list(
    longitude = mean(longitude, na.rm = TRUE),
    latitude = mean(latitude, na.rm = TRUE)
  )

  if (!is.null(location_col)) {
    row[[location_col]] <- unique(as.character(group_data[[location_col]]))[1]

    unique_pairs <- unique(data.frame(longitude = longitude, latitude = latitude))
    if (nrow(unique_pairs) > 1L) {
      warning(
        sprintf(
          "Multiple coordinate pairs detected within location '%s'; returning the centroid.",
          row[[location_col]]
        ),
        call. = FALSE
      )
    }

    row <- c(
      setNames(list(row[[location_col]]), location_col),
      list(longitude = row$longitude, latitude = row$latitude)
    )
    return(row)
  }

  row
}

#' Summarize an auxiliary metadata column at location level
#'
#' @param x Vector of values.
#'
#' @return A scalar summary value.
#' @keywords internal
summarize_extra_column <- function(x) {
  if (is.numeric(x)) {
    return(mean(x, na.rm = TRUE))
  }

  unique_values <- unique(as.character(x[!is.na(x)]))
  if (length(unique_values) == 0L) {
    return(NA_character_)
  }

  paste(unique_values, collapse = ";")
}

#' Transform a diversity trace statistic
#'
#' @param x Numeric scalar.
#' @param transform Transformation name.
#'
#' @return Numeric scalar.
#' @keywords internal
transform_div_trace <- function(x, transform) {
  switch(
    transform,
    none = x,
    log = log(x),
    log1p = log1p(x)
  )
}
