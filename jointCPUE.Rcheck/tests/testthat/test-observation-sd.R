test_that("jointCPUE can estimate fleet-specific observation SD", {
  data_input <- data.frame(
    cpue = c(1.0, 1.1, 0.9, 4.0, 7.0, 10.0, 1.2, 0.8, 1.0, 3.5, 8.5, 11.0),
    lon = c(150.0, 150.2, 150.4, 150.1, 150.3, 150.5, 151.0, 151.2, 151.4, 151.1, 151.3, 151.5),
    lat = c(40.0, 40.1, 40.2, 40.3, 40.4, 40.5, 41.0, 41.1, 41.2, 41.3, 41.4, 41.5),
    tid = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
    fleetid = c(0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1)
  )

  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  fit <- jointCPUE(
    data_utm = utm$data_utm,
    mesh = mesh,
    pop_spatial = "off",
    pop_spatiotemporal = "off",
    q_diffs_system = "off",
    q_diffs_time = "off",
    q_diffs_spatial = "off",
    index = "monthly",
    obs_sd = "fleet",
    silent = TRUE
  )

  rep_obj <- fit$obj$report(fit$obj$env$last.par.best)

  expect_equal(fit$settings$obs_sd, "fleet")
  expect_equal(fit$data_tmb$use_fleet_sd, 1L)
  expect_equal(length(rep_obj$sd_fleet), fit$data_tmb$n_f)
  expect_true(all(is.finite(rep_obj$sd_fleet)))
  expect_true(all(rep_obj$sd_fleet > 0))
})
