test_that("disabled systematic fleet effects add no objective constant", {
  data_input <- tiny_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )
  prep <- make_data(utm$data_utm, mesh, index = "monthly")

  objective_with_fixed_fleet_sd <- function(fleet_sd) {
    data_tmb <- prep$data
    data_tmb$use_pop_spatial <- 0L
    data_tmb$use_pop_spatiotemporal <- 0L
    data_tmb$use_pop_spatiotemporal_rw <- 0L
    data_tmb$use_month_fe <- 0L
    data_tmb$use_fleet_sd <- 0L
    data_tmb$use_q_diffs_system <- 0L
    data_tmb$use_q_diffs_time <- 0L
    data_tmb$use_q_diffs_spatial <- 0L

    parameters <- jointCPUE:::.make_parameters_jointCPUE(
      n_t = data_tmb$n_t,
      n_f = data_tmb$n_f,
      n_s = data_tmb$spde$n_s,
      n_m = data_tmb$n_m
    )
    parameters$fleet_ln_std_dev <- log(fleet_sd)
    map <- jointCPUE:::.make_map_jointCPUE(
      parameters = parameters,
      n_f = data_tmb$n_f,
      pop_spatial = "off",
      pop_spatiotemporal = "off",
      q_diffs_system = "off",
      q_diffs_time = "off",
      q_diffs_spatial = "off",
      index = "monthly",
      month_diffs = "off",
      obs_sd = "shared",
      n_m = data_tmb$n_m,
      has_tf = data_tmb$has_tf
    )
    obj <- TMB::MakeADFun(
      data = data_tmb,
      parameters = parameters,
      map = map,
      random = character(0),
      DLL = "jointCPUE",
      silent = TRUE
    )
    obj$fn(obj$par)
  }

  expect_equal(
    objective_with_fixed_fleet_sd(1),
    objective_with_fixed_fleet_sd(2),
    tolerance = 1e-10
  )
})
