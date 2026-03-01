tiny_cpue_data <- function(cpue = c(1, 2, 3, 4), month = c(5, 6, 7, 8)) {
  data.frame(
    cpue = cpue,
    lon = c(150, 150.5, 151, 151.5),
    lat = c(40, 40.5, 41, 41.5),
    tid = c(0, 0, 1, 1),
    fleetid = c(0, 1, 0, 1),
    month = month
  )
}

mock_jointCPUE_fit <- function(
    q_diffs_system = "on",
    q_diffs_time = "on",
    q_diffs_spatial = "on",
    index = "yearly",
    month_diffs = "on",
    month_values = c(4L, 6L, 11L)
) {
  month_effect_m <- if (length(month_values) <= 1L) {
    0
  } else {
    c(seq_len(length(month_values) - 1L) * 0.1, -sum(seq_len(length(month_values) - 1L) * 0.1))
  }

  report_out <- list(
    fleet_f = c(0.4, -0.2),
    fleet_t = matrix(c(1, 2, 3, 10, 20, 30), nrow = 3, ncol = 2),
    fleet_s = matrix(c(0.1, 0.2, -0.1, -0.2), nrow = 2, ncol = 2),
    month_effect_m = month_effect_m
  )

  obj <- list(
    env = list(last.par.best = 0),
    report = function(par) report_out
  )

  structure(
    list(
      obj = obj,
      data_tmb = list(
        n_f = 3L,
        n_t = 3L,
        has_tf = matrix(c(1L, 1L, 0L, 1L, 0L, 1L), nrow = 3, ncol = 2),
        A_gs = Matrix::Matrix(diag(2), sparse = TRUE)
      ),
      prep = list(
        key = data.frame(
          utm_x_scale = c(0, 1),
          utm_y_scale = c(0, 1),
          area_km2 = c(100, 100),
          area_km2_scaled = c(1, 1)
        )
      ),
      settings = list(
        q_diffs_system = q_diffs_system,
        q_diffs_time = q_diffs_time,
        q_diffs_spatial = q_diffs_spatial,
        index = index,
        month_diffs = month_diffs,
        month_values = month_values
      )
    ),
    class = "jointCPUE"
  )
}
