#' Calculate derived polygon metrics from axis summaries
#'
#' Combines polygon-level summaries of modeled ordination-axis surfaces into
#' interpretable reporting metrics. These values summarize predicted
#' ordination-space surfaces; they are not direct estimates of nucleotide
#' diversity.
#'
#' @param polygon_summaries Data frame returned by
#'   [aggregate_axis_rasters_to_polygons()].
#' @param axis_prefixes Optional axis prefixes or complete axis names to retain.
#' @param axis_numbers Optional integer axis numbers to retain.
#' @param pred_stat Raster statistic for predicted axis means.
#' @param heterogeneity_stat Raster statistic for within-polygon axis
#'   heterogeneity.
#' @param uncertainty_stat Raster statistic for posterior uncertainty.
#' @param polygon_id_col Polygon identifier column, when present.
#'
#' @return The input data frame with derived polygon metrics appended.
#' @export
calculate_polygon_axis_metrics <- function(
  polygon_summaries,
  axis_prefixes = NULL,
  axis_numbers = NULL,
  pred_stat = "mean",
  heterogeneity_stat = "sd",
  uncertainty_stat = "mean",
  polygon_id_col = "polygon_id"
) {
  if (!is.data.frame(polygon_summaries)) {
    stop("'polygon_summaries' must be a data frame.", call. = FALSE)
  }
  validate_polygon_metric_options(
    axis_prefixes = axis_prefixes,
    axis_numbers = axis_numbers,
    pred_stat = pred_stat,
    heterogeneity_stat = heterogeneity_stat,
    uncertainty_stat = uncertainty_stat,
    polygon_id_col = polygon_id_col
  )

  axis_table <- detect_polygon_metric_axes(
    names(polygon_summaries),
    pred_stat = pred_stat,
    axis_prefixes = axis_prefixes,
    axis_numbers = axis_numbers
  )
  if (nrow(axis_table) == 0L) {
    stop(
      "No usable predicted axis mean columns were found in 'polygon_summaries'.",
      call. = FALSE
    )
  }

  predicted_columns <- axis_table$prediction_column
  predicted_values <- as.matrix(polygon_summaries[predicted_columns])
  storage.mode(predicted_values) <- "numeric"
  axes_present <- !is.na(predicted_values)

  output <- polygon_summaries
  output$n_axes_used <- rowSums(axes_present)
  output$axis_centroid_distance <- sqrt(row_sum_of_squares(predicted_values))

  heterogeneity_columns <- paste0(
    axis_table$axis_name,
    "_mean_pred_",
    heterogeneity_stat
  )
  heterogeneity_values <- extract_polygon_metric_columns(
    polygon_summaries,
    heterogeneity_columns
  )
  output$axis_surface_dispersion <- row_sum_of_squares(heterogeneity_values)

  uncertainty_columns <- paste0(
    axis_table$axis_name,
    "_mean_sd_",
    uncertainty_stat
  )
  uncertainty_values <- extract_polygon_metric_columns(
    polygon_summaries,
    uncertainty_columns
  )
  output$mean_axis_uncertainty <- row_mean_or_na(uncertainty_values)

  output
}

#' Validate options for polygon metrics
#'
#' @param axis_prefixes Optional prefixes.
#' @param axis_numbers Optional axis numbers.
#' @param pred_stat Prediction statistic.
#' @param heterogeneity_stat Heterogeneity statistic.
#' @param uncertainty_stat Uncertainty statistic.
#' @param polygon_id_col Optional identifier column.
#'
#' @return TRUE, invisibly.
#' @keywords internal
validate_polygon_metric_options <- function(
  axis_prefixes,
  axis_numbers,
  pred_stat,
  heterogeneity_stat,
  uncertainty_stat,
  polygon_id_col
) {
  stats <- c(pred_stat, heterogeneity_stat, uncertainty_stat)
  if (!is.character(stats) || any(lengths(stats) != 1L) ||
      any(is.na(stats)) || any(!nzchar(stats))) {
    stop("Statistic arguments must be non-empty character strings.", call. = FALSE)
  }
  if (!is.null(axis_prefixes) && (!is.character(axis_prefixes) ||
      length(axis_prefixes) == 0L || any(is.na(axis_prefixes)) ||
      any(!nzchar(axis_prefixes)))) {
    stop("'axis_prefixes' must be a non-empty character vector when supplied.", call. = FALSE)
  }
  if (!is.null(axis_numbers) && (!is.numeric(axis_numbers) ||
      length(axis_numbers) == 0L || any(is.na(axis_numbers)) ||
      any(axis_numbers < 1) || any(axis_numbers != as.integer(axis_numbers)))) {
    stop("'axis_numbers' must be positive integers when supplied.", call. = FALSE)
  }
  if (!is.null(polygon_id_col) && (!is.character(polygon_id_col) ||
      length(polygon_id_col) != 1L || is.na(polygon_id_col) ||
      !nzchar(polygon_id_col))) {
    stop("'polygon_id_col' must be a non-empty character string or NULL.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Detect polygon-summary axes from prediction columns
#'
#' @param column_names Polygon summary column names.
#' @param pred_stat Prediction statistic.
#' @param axis_prefixes Optional prefixes.
#' @param axis_numbers Optional axis numbers.
#'
#' @return A data frame of detected axes.
#' @keywords internal
detect_polygon_metric_axes <- function(
  column_names,
  pred_stat,
  axis_prefixes = NULL,
  axis_numbers = NULL
) {
  pattern <- paste0("^((PC|LD)_?([0-9]+))_mean_pred_", escape_regex(pred_stat), "$")
  matches <- regexec(pattern, column_names)
  captures <- regmatches(column_names, matches)
  keep <- lengths(captures) == 4L
  if (!any(keep)) {
    return(data.frame(
      axis_name = character(0),
      axis_prefix = character(0),
      axis_number = integer(0),
      prediction_column = character(0),
      stringsAsFactors = FALSE
    ))
  }

  captures <- captures[keep]
  axis_table <- data.frame(
    axis_name = vapply(captures, function(match) match[2], character(1)),
    axis_prefix = vapply(captures, function(match) match[3], character(1)),
    axis_number = as.integer(vapply(captures, function(match) match[4], character(1))),
    prediction_column = column_names[keep],
    stringsAsFactors = FALSE
  )
  if (!is.null(axis_prefixes)) {
    axis_table <- axis_table[
      axis_table$axis_prefix %in% axis_prefixes | axis_table$axis_name %in% axis_prefixes,
      ,
      drop = FALSE
    ]
  }
  if (!is.null(axis_numbers)) {
    axis_table <- axis_table[axis_table$axis_number %in% axis_numbers, , drop = FALSE]
  }
  axis_table[order(axis_table$axis_prefix, axis_table$axis_number), , drop = FALSE]
}

#' Escape a literal string for a regular expression
#'
#' @param value Character scalar.
#'
#' @return Character scalar.
#' @keywords internal
escape_regex <- function(value) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\1", value)
}

#' Extract selected metric columns with missing columns represented by NA
#'
#' @param data Polygon summary data frame.
#' @param columns Requested columns.
#'
#' @return Numeric matrix.
#' @keywords internal
extract_polygon_metric_columns <- function(data, columns) {
  values <- lapply(columns, function(column) {
    if (!column %in% names(data)) {
      return(rep(NA_real_, nrow(data)))
    }
    if (!is.numeric(data[[column]])) {
      stop(
        sprintf("Polygon metric column '%s' must be numeric.", column),
        call. = FALSE
      )
    }
    data[[column]]
  })
  values <- do.call(cbind, values)
  if (is.null(dim(values))) {
    values <- matrix(values, ncol = 1L)
  }
  storage.mode(values) <- "numeric"
  values
}

#' Calculate a row-wise sum of squares with all-missing rows as NA
#'
#' @param values Numeric matrix.
#'
#' @return Numeric vector.
#' @keywords internal
row_sum_of_squares <- function(values) {
  has_values <- rowSums(!is.na(values)) > 0L
  output <- rowSums(values ^ 2, na.rm = TRUE)
  output[!has_values] <- NA_real_
  output
}

#' Calculate a row-wise mean with all-missing rows as NA
#'
#' @param values Numeric matrix.
#'
#' @return Numeric vector.
#' @keywords internal
row_mean_or_na <- function(values) {
  has_values <- rowSums(!is.na(values)) > 0L
  output <- rowMeans(values, na.rm = TRUE)
  output[!has_values] <- NA_real_
  output
}
