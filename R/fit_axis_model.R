#' Fit a Gaussian SPDE model for one ordination axis
#'
#' Fits a single-axis Gaussian spatial model using `INLA` and an SPDE latent
#' field. This function is intended for location-level summary responses such
#' as `PC1_mean`, `PC2_mean`, `LD1_mean`, or `LD2_mean`.
#'
#' @param points_sf Projected `sf` points returned by [prepare_spatial_points()].
#' @param response_col Name of the numeric response column to model.
#' @param mesh Either an `INLA` mesh object or the list returned by
#'   [make_inla_mesh()].
#' @param iid_effects Optional character vector of columns in `points_sf` to
#'   include as iid random effects.
#' @param family Likelihood family passed to [INLA::inla()].
#' @param compute_criteria Logical; if `TRUE`, request WAIC and DIC when
#'   available.
#' @param verbose Logical; passed to [INLA::inla()].
#'
#' @return A list with fitted model components and a compact model summary.
#' @export
fit_axis_model <- function(
  points_sf,
  response_col,
  mesh,
  iid_effects = NULL,
  family = "gaussian",
  compute_criteria = TRUE,
  verbose = FALSE
) {
  validate_projected_sf_points(points_sf, object_name = "points_sf")

  if (!response_col %in% names(points_sf)) {
    stop(sprintf("Response column '%s' was not found.", response_col), call. = FALSE)
  }

  response <- points_sf[[response_col]]
  if (!is.numeric(response)) {
    stop(sprintf("Response column '%s' must be numeric.", response_col), call. = FALSE)
  }

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required for `fit_axis_model()`.", call. = FALSE)
  }

  mesh_object <- normalize_inla_mesh(mesh)
  loc <- sf::st_coordinates(points_sf)[, 1:2, drop = FALSE]

  keep_rows <- !is.na(response)
  dropped_n <- sum(!keep_rows)
  if (dropped_n > 0L) {
    warning(
      sprintf("Dropped %d observations with missing response values in '%s'.", dropped_n, response_col),
      call. = FALSE
    )
  }

  points_model <- points_sf[keep_rows, , drop = FALSE]
  response <- response[keep_rows]
  loc <- loc[keep_rows, , drop = FALSE]

  if (length(response) < 3L) {
    stop("At least 3 non-missing observations are required to fit the axis model.", call. = FALSE)
  }

  effect_data <- build_iid_effect_data(points_model, iid_effects = iid_effects)
  spde <- build_axis_spde(mesh_object = mesh_object, response = response)
  spatial_index <- INLA::inla.spde.make.index(name = "spatial_field", n.spde = spde$n.spde)
  a_matrix <- INLA::inla.spde.make.A(mesh = mesh_object, loc = loc)

  stack_effects <- c(
    list(
      intercept = rep(1, nrow(points_model)),
      spatial_field = spatial_index$spatial_field
    ),
    effect_data$effects
  )

  stack_data <- c(
    list(response = response),
    effect_data$data
  )

  estimation_stack <- INLA::inla.stack(
    data = stack_data,
    A = c(list(1, a_matrix), effect_data$a_matrices),
    effects = stack_effects,
    tag = "estimation"
  )

  model_formula <- build_axis_formula(iid_effects = names(effect_data$effects))
  environment(model_formula) <- list2env(
    list(spde = spde),
    parent = parent.frame()
  )

  model_fit <- INLA::inla(
    formula = model_formula,
    family = family,
    data = INLA::inla.stack.data(estimation_stack),
    control.predictor = list(A = INLA::inla.stack.A(estimation_stack), compute = TRUE),
    control.compute = list(
      waic = isTRUE(compute_criteria),
      dic = isTRUE(compute_criteria)
    ),
    verbose = verbose
  )

  model_summary <- list(
    response_col = response_col,
    n_observations = nrow(points_model),
    n_mesh_vertices = nrow(mesh_object$loc[, 1:2, drop = FALSE]),
    fixed_effect_summary = model_fit$summary.fixed,
    hyperparameter_summary = model_fit$summary.hyperpar,
    waic = if (!is.null(model_fit$waic$waic)) model_fit$waic$waic else NA_real_,
    dic = if (!is.null(model_fit$dic$dic)) model_fit$dic$dic else NA_real_
  )

  list(
    model = model_fit,
    spde = spde,
    mesh = mesh_object,
    stack = estimation_stack,
    response_col = response_col,
    fitted_data = points_model,
    formula = model_formula,
    model_summary = model_summary
  )
}

#' Normalize a mesh input for axis modeling
#'
#' @param mesh A mesh object or mesh list.
#'
#' @return An INLA mesh object.
#' @keywords internal
normalize_inla_mesh <- function(mesh) {
  if (is.list(mesh) && !is.null(mesh$mesh)) {
    mesh <- mesh$mesh
  }

  if (is.null(mesh) || is.null(mesh$loc)) {
    stop("`mesh` must be an INLA mesh object or the list returned by `make_inla_mesh()`.", call. = FALSE)
  }

  mesh
}

#' Build an SPDE object for one axis model
#'
#' @param mesh_object INLA mesh object.
#' @param response Numeric response vector.
#'
#' @return An SPDE model object.
#' @keywords internal
build_axis_spde <- function(mesh_object, response) {
  vertex_loc <- mesh_object$loc[, 1:2, drop = FALSE]
  x_range <- diff(range(vertex_loc[, 1]))
  y_range <- diff(range(vertex_loc[, 2]))
  domain_scale <- max(c(x_range, y_range))
  if (!is.finite(domain_scale) || domain_scale <= 0) {
    domain_scale <- 1
  }

  response_scale <- stats::sd(response, na.rm = TRUE)
  if (!is.finite(response_scale) || response_scale <= 0) {
    response_scale <- 1
  }

  INLA::inla.spde2.pcmatern(
    mesh = mesh_object,
    alpha = 2,
    prior.range = c(domain_scale / 5, 0.5),
    prior.sigma = c(response_scale, 0.01)
  )
}

#' Build iid effect data for an axis model
#'
#' @param points_sf Projected sf points retained for modeling.
#' @param iid_effects Optional effect columns.
#'
#' @return A list with stack effects, stack data, and A matrices.
#' @keywords internal
build_iid_effect_data <- function(points_sf, iid_effects = NULL) {
  if (is.null(iid_effects) || length(iid_effects) == 0L) {
    return(list(effects = list(), data = list(), a_matrices = list()))
  }

  missing_effects <- setdiff(iid_effects, names(points_sf))
  if (length(missing_effects) > 0L) {
    stop(
      sprintf("IID effect columns were not found: %s", paste(missing_effects, collapse = ", ")),
      call. = FALSE
    )
  }

  effects <- list()
  data <- list()
  a_matrices <- list()

  for (effect_name in iid_effects) {
    effect_values <- points_sf[[effect_name]]
    if (all(is.na(effect_values))) {
      stop(sprintf("IID effect column '%s' contains only missing values.", effect_name), call. = FALSE)
    }

    effect_index <- as.integer(as.factor(effect_values))
    if (anyNA(effect_index)) {
      stop(sprintf("IID effect column '%s' contains missing values.", effect_name), call. = FALSE)
    }

    effects[[effect_name]] <- seq_len(length(unique(effect_index)))
    data[[effect_name]] <- effect_index
    a_matrices[[length(a_matrices) + 1L]] <- INLA::inla.spde.make.A(
      loc = matrix(seq_along(effect_index), ncol = 1),
      A.loc = Matrix::Diagonal(n = length(effect_index)),
      index = effect_index,
      n.spde = length(unique(effect_index))
    )
  }

  list(effects = effects, data = data, a_matrices = a_matrices)
}

#' Build the model formula for one axis fit
#'
#' @param iid_effects Character vector of iid effect names.
#'
#' @return A formula object.
#' @keywords internal
build_axis_formula <- function(iid_effects = NULL) {
  terms <- c(
    "-1",
    "intercept",
    "f(spatial_field, model = spde)"
  )

  if (!is.null(iid_effects) && length(iid_effects) > 0L) {
    terms <- c(
      terms,
      sprintf("f(%s, model = 'iid')", iid_effects)
    )
  }

  stats::as.formula(paste("response ~", paste(terms, collapse = " + ")))
}
