#' Run PCA on an encoded sequence matrix
#'
#' Performs principal components analysis on the variant matrix returned by
#' [alignment_to_variant_matrix()] and returns tidy scores linked to sequence
#' identifiers.
#'
#' @param x A numeric variant matrix, or a list containing a `variant_matrix`
#'   element.
#' @param metadata Optional data frame to join onto the score table.
#' @param id_col Name of the sequence identifier column in `metadata`.
#' @param center Logical; passed to [stats::prcomp()].
#' @param scale. Logical; passed to [stats::prcomp()].
#' @param rank. Optional maximum number of principal components to retain.
#'
#' @return A list with `scores`, `variance_summary`, and `model`.
#' @export
run_sequence_pca <- function(
  x,
  metadata = NULL,
  id_col = "sequence_id",
  center = TRUE,
  scale. = FALSE,
  rank. = NULL
) {
  variant_matrix <- extract_variant_matrix(x)
  validate_variant_matrix(variant_matrix)

  pca <- stats::prcomp(
    x = variant_matrix,
    center = center,
    scale. = scale.,
    rank. = rank.
  )

  scores <- as.data.frame(pca$x, stringsAsFactors = FALSE)
  scores <- cbind(
    data.frame(
      sequence_id = rownames(variant_matrix),
      stringsAsFactors = FALSE
    ),
    scores
  )
  scores <- join_scores_metadata(scores, metadata = metadata, id_col = id_col)

  variance_summary <- data.frame(
    component = paste0("PC", seq_along(pca$sdev)),
    standard_deviation = pca$sdev,
    variance = pca$sdev^2,
    proportion_variance = (pca$sdev^2) / sum(pca$sdev^2),
    stringsAsFactors = FALSE
  )
  variance_summary$cumulative_variance <- cumsum(variance_summary$proportion_variance)

  list(
    scores = scores,
    variance_summary = variance_summary,
    model = pca
  )
}

#' Extract a variant matrix from supported inputs
#'
#' @param x A matrix or list containing `variant_matrix`.
#'
#' @return A numeric matrix.
#' @keywords internal
extract_variant_matrix <- function(x) {
  if (is.list(x) && !is.null(x$variant_matrix)) {
    x <- x$variant_matrix
  }

  if (!is.matrix(x)) {
    stop("Expected a numeric matrix or an object with `$variant_matrix`.", call. = FALSE)
  }

  storage.mode(x) <- "double"
  x
}

#' Validate a variant matrix before ordination
#'
#' @param variant_matrix A numeric matrix.
#'
#' @return `TRUE`, invisibly.
#' @keywords internal
validate_variant_matrix <- function(variant_matrix) {
  if (!is.numeric(variant_matrix)) {
    stop("`variant_matrix` must be numeric.", call. = FALSE)
  }

  if (nrow(variant_matrix) < 1L) {
    stop("`variant_matrix` must contain at least one sequence.", call. = FALSE)
  }

  if (ncol(variant_matrix) < 1L) {
    stop(
      paste(
        "`variant_matrix` must contain at least one encoded feature.",
        "Likely causes include all sites being filtered from the alignment encoder."
      ),
      call. = FALSE
    )
  }

  if (is.null(rownames(variant_matrix)) || anyNA(rownames(variant_matrix))) {
    stop("`variant_matrix` must have sequence IDs as row names.", call. = FALSE)
  }

  if (any(!is.finite(variant_matrix))) {
    stop("`variant_matrix` must contain only finite numeric values.", call. = FALSE)
  }

  column_variances <- apply(variant_matrix, 2, stats::var)
  column_variances[is.na(column_variances)] <- 0
  if (!any(column_variances > 0)) {
    stop(
      paste(
        "`variant_matrix` must have nonzero variance in at least one column.",
        "This often indicates all retained sites encode identical sequence patterns."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Join score tables to metadata
#'
#' @param scores A data frame with a `sequence_id` column.
#' @param metadata Optional metadata data frame.
#' @param id_col Name of the sequence identifier column in `metadata`.
#'
#' @return A data frame.
#' @keywords internal
join_scores_metadata <- function(scores, metadata = NULL, id_col = "sequence_id") {
  if (is.null(metadata)) {
    return(scores)
  }

  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data frame when supplied.", call. = FALSE)
  }

  if (!id_col %in% names(metadata)) {
    stop(
      sprintf("Metadata is missing the sequence ID column '%s'.", id_col),
      call. = FALSE
    )
  }

  metadata_copy <- metadata
  names(metadata_copy)[names(metadata_copy) == id_col] <- "sequence_id"
  merge(scores, metadata_copy, by = "sequence_id", all.x = TRUE, sort = FALSE)
}
