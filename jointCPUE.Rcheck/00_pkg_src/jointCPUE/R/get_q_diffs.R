#' Extract q-difference components
#'
#' Extracts fleet-specific systematic, temporal, and spatial q-difference
#' components from a fitted `jointCPUE` model.
#'
#' Temporal deviations are centered exactly as in the C++ model:
#' each fleet-specific time series is centered by its mean over observed
#' time steps (`has_tf == 1`), and unobserved time steps are returned as `NA`.
#'
#' Spatial deviations are projected from mesh vertices onto the extrapolation
#' grid stored in the fitted object.
#'
#' @param object A fitted `jointCPUE` object.
#' @param fleet_labels Optional character vector of length `n_f` giving labels
#'   for fleet ids `0:(n_f-1)`. If `NULL`, default labels are used.
#' @param time_values Optional vector of length `n_t` used for plotting or
#'   reporting time labels. If `NULL`, `1:n_t` is used.
#'
#' @return A named list with elements `system`, `time`, and `spatial`.
#'   Components that are switched off are returned as `NULL`.
#' @export
get_q_diffs <- function(object, fleet_labels = NULL, time_values = NULL) {
  if (!inherits(object, "jointCPUE")) {
    stop("`object` must be a fitted `jointCPUE` object.", call. = FALSE)
  }

  n_f <- object$data_tmb$n_f
  n_t <- object$data_tmb$n_t

  if (is.null(fleet_labels)) {
    fleet_labels <- c("Reference", paste("Fleet", seq_len(n_f - 1L)))
  }
  if (!is.character(fleet_labels) || length(fleet_labels) != n_f) {
    stop("`fleet_labels` must be NULL or a character vector of length n_f.", call. = FALSE)
  }

  if (is.null(time_values)) {
    time_values <- seq_len(n_t)
  }
  if (length(time_values) != n_t) {
    stop("`time_values` must have length n_t.", call. = FALSE)
  }

  par_best <- object$obj$env$last.par.best
  if (is.null(par_best)) {
    stop("Could not find `last.par.best` in fitted object.", call. = FALSE)
  }
  rep_obj <- object$obj$report(par_best)

  out <- list(system = NULL, time = NULL, spatial = NULL)

  if (identical(object$settings$q_diffs_system, "on") && n_f > 1L) {
    effect <- c(0, rep_obj$fleet_f)
    out$system <- data.frame(
      fleetid = 0:(n_f - 1L),
      fleet = fleet_labels,
      effect = effect,
      reference = c(TRUE, rep(FALSE, n_f - 1L))
    )
  }

  if (identical(object$settings$q_diffs_time, "on") && n_f > 1L) {
    has_tf <- object$data_tmb$has_tf > 0L
    fleet_t <- rep_obj$fleet_t
    centered <- fleet_t

    for (j in seq_len(ncol(fleet_t))) {
      keep <- has_tf[, j]
      if (any(keep)) {
        centered[, j] <- fleet_t[, j] - mean(fleet_t[keep, j])
      } else {
        centered[, j] <- NA_real_
      }
      centered[!keep, j] <- NA_real_
    }

    time_list <- vector("list", n_f - 1L)
    for (j in seq_len(n_f - 1L)) {
      time_list[[j]] <- data.frame(
        tid = 0:(n_t - 1L),
        time = time_values,
        fleetid = j,
        fleet = fleet_labels[j + 1L],
        effect = centered[, j],
        observed = has_tf[, j]
      )
    }
    out$time <- do.call(rbind, time_list)
    rownames(out$time) <- NULL
  }

  if (identical(object$settings$q_diffs_spatial, "on") && n_f > 1L) {
    key <- object$prep$key
    proj <- as.matrix(object$data_tmb$A_gs %*% rep_obj$fleet_s)

    spatial_list <- vector("list", n_f - 1L)
    for (j in seq_len(n_f - 1L)) {
      spatial_list[[j]] <- data.frame(
        grid_id = seq_len(nrow(key)),
        fleetid = j,
        fleet = fleet_labels[j + 1L],
        utm_x_scale = key$utm_x_scale,
        utm_y_scale = key$utm_y_scale,
        area_km2 = key$area_km2,
        area_km2_scaled = key$area_km2_scaled,
        effect = proj[, j]
      )
    }
    out$spatial <- do.call(rbind, spatial_list)
    rownames(out$spatial) <- NULL
  }

  out
}
