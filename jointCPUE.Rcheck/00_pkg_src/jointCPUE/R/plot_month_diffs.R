#' Plot month fixed effects
#'
#' Plots estimated calendar-month fixed effects from a fitted `jointCPUE`
#' model. The x-axis shows only observed calendar months from the fitted data.
#'
#' @param object A fitted `jointCPUE` object.
#'
#' @return A ggplot object.
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
