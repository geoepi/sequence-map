#' Read an alignment from a FASTA file
#'
#' Reads a FASTA alignment into memory and verifies that sequence identifiers
#' are available for downstream joins to metadata.
#'
#' @param path Path to a FASTA file.
#'
#' @return A named list of sequences as returned by [ape::read.FASTA()].
#' @export
read_alignment <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Alignment file does not exist: %s", path), call. = FALSE)
  }

  alignment <- tryCatch(
    ape::read.FASTA(file = path),
    error = function(error) {
      stop(
        sprintf("Failed to read alignment from '%s': %s", path, error$message),
        call. = FALSE
      )
    }
  )

  sequence_ids <- names(alignment)
  if (length(alignment) == 0L || is.null(sequence_ids) || anyNA(sequence_ids)) {
    stop(
      sprintf("Alignment '%s' did not contain usable sequence IDs.", path),
      call. = FALSE
    )
  }

  alignment
}

#' Read sequence metadata
#'
#' Reads metadata from a delimited text file or RDS file.
#'
#' @param path Path to a metadata file. Supported extensions are `.csv`,
#'   `.tsv`, `.txt`, and `.rds`.
#'
#' @return A data frame containing sequence metadata.
#' @export
read_metadata <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Metadata file does not exist: %s", path), call. = FALSE)
  }

  extension <- tolower(tools::file_ext(path))
  metadata <- switch(
    extension,
    csv = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    tsv = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    txt = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    rds = readRDS(path),
    stop(
      sprintf(
        "Unsupported metadata format for '%s'. Use .csv, .tsv, .txt, or .rds.",
        path
      ),
      call. = FALSE
    )
  )

  if (!is.data.frame(metadata)) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  }

  metadata
}

#' Read a spatial boundary
#'
#' Reads a boundary dataset with [sf::st_read()] and confirms that the result
#' is an `sf` object.
#'
#' @param path Path to a boundary dataset readable by GDAL.
#' @param layer Optional layer name for multi-layer data sources.
#' @param quiet Logical; passed to [sf::st_read()].
#'
#' @return An `sf` object.
#' @export
read_boundary <- function(path, layer = NULL, quiet = TRUE) {
  if (!file.exists(path)) {
    stop(sprintf("Boundary file does not exist: %s", path), call. = FALSE)
  }

  boundary <- tryCatch(
    if (is.null(layer)) {
      sf::st_read(path, quiet = quiet)
    } else {
      sf::st_read(dsn = path, layer = layer, quiet = quiet)
    },
    error = function(error) {
      stop(
        sprintf("Failed to read boundary as sf from '%s': %s", path, error$message),
        call. = FALSE
      )
    }
  )

  if (!inherits(boundary, "sf")) {
    stop(
      sprintf("Boundary '%s' was read, but the result is not an sf object.", path),
      call. = FALSE
    )
  }

  boundary
}
