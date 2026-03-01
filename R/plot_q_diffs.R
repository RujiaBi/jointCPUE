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

#' Plot systematic q-differences
#'
#' Plots fleet-specific systematic q-differences relative to the reference fleet.
#'
#' @param object A fitted `jointCPUE` object.
#' @param fleet_labels Optional character vector of length `n_f`.
#'
#' @return A ggplot object.
#' @export
plot_q_diffs_system <- function(object, fleet_labels = NULL) {
  qd <- get_q_diffs(object, fleet_labels = fleet_labels)$system
  if (is.null(qd)) {
    stop("Systematic q-differences are not available. Fit with `q_diffs_system = 'on'`.", call. = FALSE)
  }

  qd$fleet <- factor(qd$fleet, levels = qd$fleet)

  ggplot2::ggplot(qd, ggplot2::aes(x = fleet, y = effect, colour = reference)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.6) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = fleet, y = 0, yend = effect),
      linewidth = 1.0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#B45309", "FALSE" = "#0F766E"), guide = "none") +
    ggplot2::labs(x = NULL, y = "Systematic q-difference") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
}

#' Plot temporal q-differences
#'
#' Plots fleet-specific temporal q-differences over time. Effects are centered
#' as in the fitted model and shown only at observed fleet-time combinations.
#'
#' @param object A fitted `jointCPUE` object.
#' @param fleet_labels Optional character vector of length `n_f`.
#' @param time_values Optional vector of length `n_t` used as the x-axis.
#'
#' @return A ggplot object.
#' @export
plot_q_diffs_time <- function(object, fleet_labels = NULL, time_values = NULL) {
  qd <- get_q_diffs(object, fleet_labels = fleet_labels, time_values = time_values)$time
  if (is.null(qd)) {
    stop("Temporal q-differences are not available. Fit with `q_diffs_time = 'on'`.", call. = FALSE)
  }

  qd_obs <- qd[qd$observed, , drop = FALSE]

  ggplot2::ggplot(qd_obs, ggplot2::aes(x = time, y = effect)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.6) +
    ggplot2::geom_line(colour = "#1D4ED8", linewidth = 0.9) +
    ggplot2::geom_point(colour = "#1D4ED8", fill = "white", shape = 21, stroke = 0.6, size = 1.8) +
    ggplot2::facet_wrap(~ fleet, scales = "free_y") +
    ggplot2::labs(x = "Time", y = "Temporal q-difference") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80")
    )
}

#' Plot spatial q-differences
#'
#' Plots fleet-specific spatial q-differences projected onto the extrapolation
#' grid used for index calculation.
#'
#' @param object A fitted `jointCPUE` object.
#' @param fleet_labels Optional character vector of length `n_f`.
#' @param point_size Point size used for grid-cell plotting.
#'
#' @return A ggplot object.
#' @export
plot_q_diffs_spatial <- function(object, fleet_labels = NULL, point_size = 2.5) {
  qd <- get_q_diffs(object, fleet_labels = fleet_labels)$spatial
  if (is.null(qd)) {
    stop("Spatial q-differences are not available. Fit with `q_diffs_spatial = 'on'`.", call. = FALSE)
  }

  ggplot2::ggplot(
    qd,
    ggplot2::aes(x = utm_x_scale, y = utm_y_scale, colour = effect)
  ) +
    ggplot2::geom_point(shape = 15, size = point_size) +
    ggplot2::facet_wrap(~ fleet) +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      name = "Spatial q-difference"
    ) +
    ggplot2::labs(x = "UTM X", y = "UTM Y") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80")
    )
}
