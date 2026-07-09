#' Run DAPC on an encoded sequence matrix
#'
#' Runs discriminant analysis of principal components using `adegenet`.
#' Groups must come from a metadata column or be inferred with
#' [adegenet::find.clusters()].
#'
#' @param x A numeric variant matrix, or a list containing a `variant_matrix`
#'   element.
#' @param metadata Optional data frame containing sequence-level metadata.
#' @param id_col Name of the sequence identifier column in `metadata`.
#' @param grouping_col Optional metadata column defining groups for DAPC.
#' @param infer_groups Logical; if `TRUE`, infer groups with
#'   [adegenet::find.clusters()] when `grouping_col` is not supplied.
#' @param n_pca Number of principal components retained before DAPC.
#' @param n_da Number of discriminant axes to retain.
#' @param n_clust Optional number of inferred clusters when
#'   `infer_groups = TRUE`.
#' @param max_n_clust Maximum clusters considered when inferring groups.
#' @param center Logical; passed to adegenet methods.
#' @param scale Logical; passed to adegenet methods.
#' @param seed Optional random seed for inferred grouping.
#'
#' @return A list with `scores`, `eigen_summary`, `group_summary`, and `model`.
#' @export
run_sequence_dapc <- function(
  x,
  metadata = NULL,
  id_col = "sequence_id",
  grouping_col = NULL,
  infer_groups = FALSE,
  n_pca = NULL,
  n_da = NULL,
  n_clust = NULL,
  max_n_clust = NULL,
  center = TRUE,
  scale = FALSE,
  seed = NULL
) {
  variant_matrix <- extract_variant_matrix(x)
  validate_variant_matrix(variant_matrix)

  grouping_result <- resolve_dapc_groups(
    variant_matrix = variant_matrix,
    metadata = metadata,
    id_col = id_col,
    grouping_col = grouping_col,
    infer_groups = infer_groups,
    n_pca = n_pca,
    n_clust = n_clust,
    max_n_clust = max_n_clust,
    center = center,
    scale = scale,
    seed = seed
  )

  if (is.null(n_pca)) {
    n_pca <- min(ncol(variant_matrix), max(1L, nrow(variant_matrix) - 1L))
  }

  group_count <- nlevels(grouping_result$groups)
  max_da <- max(1L, min(group_count - 1L, nrow(variant_matrix) - 1L))
  if (is.null(n_da)) {
    n_da <- max_da
  } else {
    n_da <- min(n_da, max_da)
  }

  dapc_fit <- adegenet::dapc(
    x = as.data.frame(variant_matrix),
    grp = grouping_result$groups,
    n.pca = n_pca,
    n.da = n_da,
    center = center,
    scale = scale
  )

  score_table <- as.data.frame(dapc_fit$ind.coord, stringsAsFactors = FALSE)
  if (ncol(score_table) == 0L) {
    score_table <- data.frame(LD1 = rep(0, nrow(variant_matrix)))
  }

  names(score_table) <- paste0("LD", seq_len(ncol(score_table)))
  scores <- cbind(
    data.frame(
      sequence_id = rownames(variant_matrix),
      group = as.character(grouping_result$groups),
      assigned_group = as.character(dapc_fit$assign),
      stringsAsFactors = FALSE
    ),
    score_table
  )
  scores <- join_scores_metadata(scores, metadata = metadata, id_col = id_col)

  posterior <- as.data.frame(dapc_fit$posterior, stringsAsFactors = FALSE)
  names(posterior) <- paste0("posterior_", names(posterior))
  scores <- cbind(scores, posterior)

  eigen_summary <- data.frame(
    discriminant = paste0("LD", seq_along(dapc_fit$eig)),
    eigenvalue = dapc_fit$eig,
    proportion_eigenvalue = dapc_fit$eig / sum(dapc_fit$eig),
    stringsAsFactors = FALSE
  )
  eigen_summary$cumulative_eigenvalue <- cumsum(eigen_summary$proportion_eigenvalue)

  group_summary <- data.frame(
    group = names(table(grouping_result$groups)),
    n_sequences = as.integer(table(grouping_result$groups)),
    stringsAsFactors = FALSE
  )

  list(
    scores = scores,
    eigen_summary = eigen_summary,
    group_summary = group_summary,
    inferred_groups = grouping_result$inferred,
    model = dapc_fit
  )
}

#' Resolve DAPC grouping assignments
#'
#' @param variant_matrix Numeric variant matrix.
#' @param metadata Optional metadata data frame.
#' @param id_col Name of the sequence identifier column in `metadata`.
#' @param grouping_col Optional metadata grouping column.
#' @param infer_groups Logical; whether to infer groups.
#' @param n_pca Number of PCs for clustering.
#' @param n_clust Optional fixed cluster count.
#' @param max_n_clust Maximum clusters considered.
#' @param center Logical; passed to adegenet.
#' @param scale Logical; passed to adegenet.
#' @param seed Optional random seed.
#'
#' @return A list with `groups` and `inferred`.
#' @keywords internal
resolve_dapc_groups <- function(
  variant_matrix,
  metadata,
  id_col,
  grouping_col,
  infer_groups,
  n_pca,
  n_clust,
  max_n_clust,
  center,
  scale,
  seed
) {
  if (!is.null(grouping_col)) {
    if (is.null(metadata)) {
      stop("`metadata` must be supplied when `grouping_col` is used.", call. = FALSE)
    }

    if (!is.data.frame(metadata)) {
      stop("`metadata` must be a data frame.", call. = FALSE)
    }

    required_columns <- c(id_col, grouping_col)
    missing_columns <- setdiff(required_columns, names(metadata))
    if (length(missing_columns) > 0L) {
      stop(
        sprintf(
          "Metadata is missing required columns for DAPC: %s",
          paste(missing_columns, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    match_index <- match(rownames(variant_matrix), metadata[[id_col]])
    missing_ids <- rownames(variant_matrix)[is.na(match_index)]
    if (length(missing_ids) > 0L) {
      stop(
        sprintf(
          "Sequence IDs in the variant matrix are missing from metadata (%d missing). Examples: %s",
          length(missing_ids),
          format_example_ids(missing_ids)
        ),
        call. = FALSE
      )
    }

    groups <- metadata[[grouping_col]][match_index]
    if (anyNA(groups) || any(!nzchar(as.character(groups)))) {
      bad_ids <- rownames(variant_matrix)[is.na(groups) | !nzchar(as.character(groups))]
      stop(
        sprintf(
          "Grouping assignments are missing for some sequences. Examples: %s",
          format_example_ids(bad_ids)
        ),
        call. = FALSE
      )
    }

    groups <- as.factor(groups)
    if (nlevels(groups) < 2L) {
      stop("DAPC requires at least two groups.", call. = FALSE)
    }

    return(list(groups = groups, inferred = FALSE))
  }

  if (!isTRUE(infer_groups)) {
    stop(
      "Provide `grouping_col` or set `infer_groups = TRUE` to run DAPC.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(n_pca)) {
    n_pca <- min(ncol(variant_matrix), max(1L, nrow(variant_matrix) - 1L))
  }

  if (is.null(max_n_clust)) {
    max_n_clust <- max(2L, min(10L, nrow(variant_matrix) - 1L))
  }

  clustering <- adegenet::find.clusters(
    x = as.data.frame(variant_matrix),
    n.pca = n_pca,
    n.clust = n_clust,
    choose.n.clust = FALSE,
    criterion = "min",
    max.n.clust = max_n_clust,
    center = center,
    scale = scale
  )

  groups <- as.factor(clustering$grp)
  if (nlevels(groups) < 2L) {
    stop("Inferred grouping produced fewer than two groups for DAPC.", call. = FALSE)
  }

  list(groups = groups, inferred = TRUE)
}
