#' sequencemap: Spatial Mapping of Viral Sequence Variation
#'
#' Tools for converting aligned viral sequence data and sampling coordinates
#' into PCA-first ordination summaries, optional INLA-SPDE spatial predictions,
#' raster outputs, and polygon reporting summaries.
#'
#' ## Input contract
#'
#' * FASTA alignments must contain aligned, same-length DNA sequences with
#'   unique sequence IDs. Gaps and ambiguous DNA states are allowed and retained
#'   as quality-control states.
#' * Metadata must contain matching sequence IDs and numeric longitude and
#'   latitude columns. Location and group columns are optional.
#' * Boundaries must be valid polygon data readable by sf, have a usable CRS,
#'   and cover sampled locations.
#'
#' Use [run_sequence_map_workflow()] for the end-to-end PCA-first workflow.
#' Use [summarize_workflow_status()] and [write_workflow_status()] to inspect
#' which optional stages succeeded, failed, or were skipped.
#'
#' @name sequencemap
#' @aliases sequencemap-package
#' @docType package
#' @keywords package
"_PACKAGE"
