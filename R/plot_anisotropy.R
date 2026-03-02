#' Plot anisotropy ellipses for the fitted spatial field
#'
#' Visualizes the anisotropy implied by `ln_H_input` together with the fitted
#' range parameter for the population spatial field.
#'
#' @param object An object of class `jointCPUE` returned by [jointCPUE()], or a
#'   TMB object list containing `obj$env$last.par.best`.
#' @param colors Character vector of length 1 or 2 giving ellipse colors. If a
#'   single color is supplied, it is recycled.
#' @param labels Character vector of length 1 or 2 used in the legend. If a
#'   single label is supplied, it is recycled.
#' @param n_points Number of points used to draw each ellipse.
#'
#' @return A ggplot object.
#' @export
plot_anisotropy <- function(
    object,
    colors = c("#2E7D32", "#111111"),
    labels = c("Population field", "Population field"),
    n_points = 361L
) {
  par_best <- NULL

  if (inherits(object, "jointCPUE")) {
    par_best <- object$obj$env$last.par.best
  } else if (is.list(object) &&
             !is.null(object$obj) &&
             !is.null(object$obj$env) &&
             !is.null(object$obj$env$last.par.best)) {
    par_best <- object$obj$env$last.par.best
  }

  if (is.null(par_best)) {
    stop(
      "`object` must be a `jointCPUE` fit or contain `obj$env$last.par.best`.",
      call. = FALSE
    )
  }

  ln_H_input <- as.numeric(par_best[grep("^ln_H_input", names(par_best))])
  if (length(ln_H_input) != 2L) {
    stop("Could not find the two anisotropy parameters `ln_H_input`.", call. = FALSE)
  }

  H <- matrix(
    c(
      exp(ln_H_input[1]), ln_H_input[2],
      ln_H_input[2], (1 + ln_H_input[2]^2) / exp(ln_H_input[1])
    ),
    nrow = 2,
    byrow = TRUE
  )

  eig <- eigen(H)

  extract_range <- function(idx) {
    nm_range <- paste0("ln_range_", idx)
    nm_kappa <- paste0("ln_kappa_", idx)

    if (nm_range %in% names(par_best)) {
      return(exp(unname(par_best[[nm_range]])))
    }
    if (nm_kappa %in% names(par_best)) {
      return(sqrt(8.0) / exp(unname(par_best[[nm_kappa]])))
    }

    stop(
      "Could not find `", nm_range, "` or `", nm_kappa, "` in `last.par.best`.",
      call. = FALSE
    )
  }

  ranges <- extract_range(1)

  if (length(colors) < 1L) {
    stop("`colors` must have length >= 1.", call. = FALSE)
  }
  if (length(labels) < 1L) {
    stop("`labels` must have length >= 1.", call. = FALSE)
  }
  colors <- rep(colors, length.out = 1L)
  labels <- rep(labels, length.out = 1L)

  n_points <- as.integer(n_points)
  if (is.na(n_points) || n_points < 20L) {
    stop("`n_points` must be an integer >= 20.", call. = FALSE)
  }

  major <- eig$vectors[, 1] * eig$values[1] * ranges
  minor <- eig$vectors[, 2] * eig$values[2] * ranges

  theta <- seq(0, 2 * pi, length.out = n_points)
  xy <- vapply(
    theta,
    function(th) major * cos(th) + minor * sin(th),
    numeric(2)
  )

  ellipse_df <- data.frame(
    x = xy[1, ],
    y = xy[2, ],
    component = labels[1]
  )

  axis_df <- data.frame(
    x = c(-major[1], major[1], -minor[1], minor[1]),
    y = c(-major[2], major[2], -minor[2], minor[2]),
    xend = c(major[1], -major[1], minor[1], -minor[1]),
    yend = c(major[2], -major[2], minor[2], -minor[2]),
    component = labels[1],
    axis_type = rep(c("major", "major", "minor", "minor"), each = 1L)
  )

  lim <- max(abs(c(ellipse_df$x, ellipse_df$y)), na.rm = TRUE) * 1.1

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = ellipse_df,
      ggplot2::aes(x = x, y = y, colour = component),
      linewidth = 1
    ) +
    ggplot2::geom_segment(
      data = axis_df,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = component),
      linewidth = 0.4,
      alpha = 0.7
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.4) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::scale_colour_manual(values = stats::setNames(colors, labels)) +
    ggplot2::labs(
      x = "Scaled Easting",
      y = "Scaled Northing",
      colour = NULL,
      title = "Distance at 10% correlation"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    )
}
