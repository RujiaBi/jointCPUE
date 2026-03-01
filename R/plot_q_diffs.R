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

  ggplot2::ggplot(qd, ggplot2::aes(x = time, y = effect)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.6) +
    ggplot2::geom_line(colour = "#1D4ED8", linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_point(colour = "#1D4ED8", fill = "white", shape = 21, stroke = 0.6, size = 1.8, na.rm = TRUE) +
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
