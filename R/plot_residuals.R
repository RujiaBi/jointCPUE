#' Plot residual diagnostics
#'
#' Creates two residual diagnostic plots:
#' \itemize{
#'   \item an observed-versus-predicted scatter plot with a 1:1 reference line
#'   \item a spatial plot of mean residuals aggregated within 2D bins
#' }
#'
#' The spatial plot uses binning by default to reduce overplotting when many
#' observations overlap in space.
#'
#' @param pred_df A data.frame returned by [get_predicted()].
#' @param observed_col Name of the observed response column. Default `"observed"`.
#' @param x_col Name of the x-coordinate column. Default `"utm_x_scale"`.
#' @param y_col Name of the y-coordinate column. Default `"utm_y_scale"`.
#' @param scatter_scale Scale for the observed-versus-predicted plot:
#'   `"log"` (default) or `"response"`.
#' @param residual Which residual to plot spatially: `"log"` or `"raw"`.
#' @param bins Number of bins for the spatial residual plot.
#'
#' @return A named list of ggplot objects with elements `observed_predicted`
#'   and `spatial_residual`.
#' @export
plot_residuals <- function(
    pred_df,
    observed_col = "observed",
    x_col = "utm_x_scale",
    y_col = "utm_y_scale",
    scatter_scale = c("log", "response"),
    residual = c("log", "raw"),
    bins = 40
) {
  scatter_scale <- match.arg(scatter_scale)
  residual <- match.arg(residual)

  req <- c(
    observed_col, x_col, y_col, "fitted", "observed_log",
    "predicted_log", "residual_raw", "residual_log"
  )
  miss <- setdiff(req, names(pred_df))
  if (length(miss)) {
    stop(
      "`pred_df` is missing required columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  z_col <- if (residual == "log") "residual_log" else "residual_raw"
  z_lab <- if (residual == "log") {
    "Mean log residual (observed - predicted)"
  } else {
    "Mean raw residual (observed - predicted)"
  }

  if (scatter_scale == "log") {
    plot_dat <- data.frame(
      observed = pred_df[["observed_log"]],
      predicted = pred_df[["predicted_log"]],
      x = pred_df[[x_col]],
      y = pred_df[[y_col]],
      z = pred_df[[z_col]]
    )
    x_lab <- "Observed log CPUE"
    y_lab <- "Predicted log CPUE"
  } else {
    plot_dat <- data.frame(
      observed = pred_df[[observed_col]],
      predicted = pred_df[["fitted"]],
      x = pred_df[[x_col]],
      y = pred_df[[y_col]],
      z = pred_df[[z_col]]
    )
    x_lab <- "Observed CPUE"
    y_lab <- "Predicted CPUE"
  }

  p_obs_pred <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(x = observed, y = predicted)
  ) +
    ggplot2::geom_point(
      colour = "#0F766E",
      alpha = 0.28,
      size = 1.6
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      colour = "#C0392B",
      linewidth = 0.9
    ) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92")
    )

  p_spatial <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(x = x, y = y, z = z)
  ) +
    ggplot2::stat_summary_2d(
      ggplot2::aes(fill = ggplot2::after_stat(value)),
      fun = mean,
      bins = bins
    ) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      name = z_lab
    ) +
    ggplot2::labs(
      x = x_col,
      y = y_col
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank()
    )

  list(
    observed_predicted = p_obs_pred,
    spatial_residual = p_spatial
  )
}
