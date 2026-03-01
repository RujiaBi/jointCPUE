#' Extract month fixed effects
#'
#' Extracts estimated calendar-month fixed effects from a fitted `jointCPUE`
#' model. Effects are returned on the log scale used in the observation model
#' and follow the same sum-to-zero constraint as the underlying TMB template.
#'
#' @param object A fitted `jointCPUE` object.
#'
#' @return For `get_month_diffs()`, a data.frame with columns `month` and
#'   `effect`, or `NULL` when month fixed effects are not available for the
#'   fitted model. For `plot_month_diffs()`, a ggplot object.
#' @export
get_month_diffs <- function(object) {
  if (!inherits(object, "jointCPUE")) {
    stop("`object` must be a fitted `jointCPUE` object.", call. = FALSE)
  }

  index <- object$settings$index
  if (is.null(index)) {
    index <- object$prep$time$index
  }

  month_values <- object$settings$month_values
  if (is.null(month_values)) {
    month_values <- object$prep$time$month_values
  }

  month_diffs <- object$settings$month_diffs
  if (is.null(month_diffs)) {
    month_diffs <- if (identical(index, "yearly") && length(month_values) > 0L) "on" else "off"
  }

  if (!identical(index, "yearly") || !identical(month_diffs, "on")) {
    return(NULL)
  }

  if (is.null(month_values) || length(month_values) < 1L) {
    stop("Could not find observed month values in fitted object.", call. = FALSE)
  }

  par_best <- object$obj$env$last.par.best
  if (is.null(par_best)) {
    stop("Could not find `last.par.best` in fitted object.", call. = FALSE)
  }
  rep_obj <- object$obj$report(par_best)

  effect <- rep_obj$month_effect_m
  if (is.null(effect)) {
    stop("Could not find reported month effects in fitted object.", call. = FALSE)
  }
  if (length(effect) != length(month_values)) {
    stop("Reported month effects do not match stored observed month values.", call. = FALSE)
  }

  out <- data.frame(
    month = as.integer(month_values),
    effect = as.numeric(effect)
  )
  rownames(out) <- NULL
  out
}

#' Plot month fixed effects
#'
#' Plots estimated calendar-month fixed effects from a fitted `jointCPUE`
#' model. The x-axis shows only observed calendar months from the fitted data.
#'
#' @param object A fitted `jointCPUE` object.
#'
#' @return For `get_month_diffs()`, a data.frame with columns `month` and
#'   `effect`, or `NULL` when month fixed effects are not available for the
#'   fitted model. For `plot_month_diffs()`, a ggplot object.
#' @export
plot_month_diffs <- function(object) {
  md <- get_month_diffs(object)
  if (is.null(md)) {
    stop(
      "Month fixed effects are not available. Fit with `index = 'yearly'` and `month_diffs = 'on'`.",
      call. = FALSE
    )
  }

  x_month <- factor(md$month, levels = md$month)
  y_effect <- md$effect

  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.6) +
    ggplot2::geom_line(
      ggplot2::aes(x = x_month, y = y_effect, group = 1),
      colour = "#B45309",
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = x_month, y = y_effect),
      colour = "#B45309",
      fill = "white",
      shape = 21,
      stroke = 0.7,
      size = 2.4
    ) +
    ggplot2::labs(x = "Month", y = "Month fixed effect") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )
}
