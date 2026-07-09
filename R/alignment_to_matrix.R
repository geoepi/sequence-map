#' Convert a DNAbin alignment to a variant matrix
#'
#' Converts a DNAbin alignment into a numeric dummy-variable matrix suitable
#' for multivariate workflows such as PCA or DAPC. Canonical states (`A`, `C`,
#' `G`, `T`) are encoded as one-hot variables for each retained site. Gap,
#' ambiguous, and missing states are recorded in QC summaries and encoded as
#' zeros across the dummy variables for the site when `missing_strategy =
#' "zero"`.
#'
#' @param alignment A DNAbin alignment with sequences in rows and sites in
#'   columns.
#' @param drop_invariant Logical; if `TRUE`, retain only sites with more than
#'   one observed canonical nucleotide across sequences.
#' @param missing_strategy Strategy for handling noncanonical states. Currently
#'   only `"zero"` is supported.
#'
#' @return A list with `variant_matrix`, `site_summary`, and
#'   `sequence_summary`.
#' @export
alignment_to_variant_matrix <- function(
  alignment,
  drop_invariant = TRUE,
  missing_strategy = c("zero")
) {
  missing_strategy <- match.arg(missing_strategy)

  if (!inherits(alignment, "DNAbin")) {
    stop("`alignment` must inherit from class 'DNAbin'.", call. = FALSE)
  }

  alignment_matrix <- as.matrix(alignment)
  if (is.null(dim(alignment_matrix))) {
    alignment_matrix <- matrix(alignment_matrix, nrow = 1L)
  }

  if (nrow(alignment_matrix) == 0L || ncol(alignment_matrix) == 0L) {
    stop("`alignment` must contain at least one sequence and one site.", call. = FALSE)
  }

  if (is.null(rownames(alignment_matrix))) {
    stop("`alignment` must have sequence IDs as row names.", call. = FALSE)
  }

  character_alignment <- ape::as.character.DNAbin(alignment)
  if (is.list(character_alignment)) {
    character_alignment <- do.call(rbind, character_alignment)
  }
  character_alignment <- toupper(as.matrix(character_alignment))
  dimnames(character_alignment) <- dimnames(alignment_matrix)

  canonical_states <- c("A", "C", "G", "T")
  gap_states <- c("-", ".")
  missing_states <- c("N", "?", NA_character_)

  make_state_mask <- function(values, states) {
    matrix(
      values %in% states,
      nrow = nrow(values),
      ncol = ncol(values),
      dimnames = dimnames(values)
    )
  }

  is_canonical <- make_state_mask(character_alignment, canonical_states)
  is_gap <- make_state_mask(character_alignment, gap_states)
  is_missing <- is.na(character_alignment) | make_state_mask(character_alignment, missing_states)
  is_ambiguous <- !(is_canonical | is_gap | is_missing)

  canonical_counts <- vapply(
    canonical_states,
    function(state) colSums(character_alignment == state, na.rm = TRUE),
    numeric(ncol(character_alignment))
  )
  canonical_counts <- t(canonical_counts)
  rownames(canonical_counts) <- canonical_states
  colnames(canonical_counts) <- seq_len(ncol(character_alignment))

  n_A <- canonical_counts["A", ]
  n_C <- canonical_counts["C", ]
  n_G <- canonical_counts["G", ]
  n_T <- canonical_counts["T", ]
  n_gap <- colSums(is_gap, na.rm = TRUE)
  n_ambiguous <- colSums(is_ambiguous, na.rm = TRUE)
  n_missing <- colSums(is_missing, na.rm = TRUE)

  observed_canonical_n <- colSums(canonical_counts > 0L)
  has_canonical <- observed_canonical_n > 0L
  is_invariant <- observed_canonical_n <= 1L
  has_noncanonical <- (n_gap + n_ambiguous + n_missing) > 0L

  site_class <- ifelse(
    observed_canonical_n >= 2L,
    "variable_canonical",
    ifelse(
      observed_canonical_n == 1L & has_noncanonical,
      "invariant_canonical_with_noncanonical",
      ifelse(
        observed_canonical_n == 1L,
        "invariant_canonical",
        ifelse(
          n_gap == nrow(character_alignment),
          "all_gap",
          ifelse(
            n_missing == nrow(character_alignment),
            "all_missing",
            ifelse(
              n_ambiguous == nrow(character_alignment),
              "all_ambiguous",
              "noncanonical_only"
            )
          )
        )
      )
    )
  )

  retained <- has_canonical & (!isTRUE(drop_invariant) | !is_invariant)

  site_index <- seq_len(ncol(character_alignment))
  site_summary <- data.frame(
    site = site_index,
    site_id = sprintf("site_%s", site_index),
    n_A = as.integer(n_A),
    n_C = as.integer(n_C),
    n_G = as.integer(n_G),
    n_T = as.integer(n_T),
    n_gap = as.integer(n_gap),
    n_ambiguous = as.integer(n_ambiguous),
    n_missing = as.integer(n_missing),
    n_observed_canonical = as.integer(observed_canonical_n),
    site_class = site_class,
    retained = retained,
    stringsAsFactors = FALSE
  )

  sequence_summary <- data.frame(
    sequence_id = rownames(character_alignment),
    n_sites = ncol(character_alignment),
    n_canonical = rowSums(is_canonical, na.rm = TRUE),
    n_gap = rowSums(is_gap, na.rm = TRUE),
    n_ambiguous = rowSums(is_ambiguous, na.rm = TRUE),
    n_missing = rowSums(is_missing, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  retained_sites <- which(retained)
  n_retained <- length(retained_sites)
  if (n_retained == 0L) {
    stop(
      paste(
        "No retained sites remain after filtering.",
        "Likely causes include all sites being invariant across canonical A/C/G/T states,",
        "or all sites containing only gaps, ambiguous states, or missing values."
      ),
      call. = FALSE
    )
  }

  n_sequences <- nrow(character_alignment)
  n_columns <- n_retained * length(canonical_states)
  variant_matrix <- matrix(
    0,
    nrow = n_sequences,
    ncol = n_columns,
    dimnames = list(rownames(character_alignment), character(0))
  )

  if (missing_strategy != "zero") {
    stop(
      sprintf("`missing_strategy = '%s'` is not implemented yet.", missing_strategy),
      call. = FALSE
    )
  }

  column_names <- character(n_columns)
  column_offset <- 1L

  for (site_position in retained_sites) {
    site_label <- sprintf("site_%s", site_position)
    site_values <- character_alignment[, site_position]

    for (state in canonical_states) {
      variant_matrix[, column_offset] <- as.numeric(site_values == state)
      column_names[column_offset] <- sprintf("%s_%s", site_label, state)
      column_offset <- column_offset + 1L
    }
  }

  colnames(variant_matrix) <- column_names

  list(
    variant_matrix = variant_matrix,
    site_summary = site_summary,
    sequence_summary = sequence_summary
  )
}
