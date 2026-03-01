test_that("core user-facing functions are exported", {
  exports <- getNamespaceExports("jointCPUE")

  expect_true(all(c(
    "jointCPUE",
    "make_utm",
    "make_mesh",
    "make_data",
    "get_index",
    "get_predicted"
  ) %in% exports))
})

test_that("jointCPUE defaults to shared observation SD", {
  data_input <- tiny_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
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
    silent = TRUE
  )

  rep_obj <- fit$obj$report(fit$obj$env$last.par.best)

  expect_equal(fit$settings$obs_sd, "shared")
  expect_equal(fit$data_tmb$use_fleet_sd, 0L)
  expect_equal(length(rep_obj$sd_fleet), fit$data_tmb$n_f)
  expect_true(all(is.finite(rep_obj$sd_fleet)))
})

test_that("jointCPUE accepts a supplied extrapolation grid", {
  data_input <- tiny_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )
  grid_common <- unique(utm$data_utm[, c("lon", "lat", "utm_x_scale", "utm_y_scale")])

  fit <- jointCPUE(
    data_utm = utm$data_utm,
    mesh = mesh,
    pop_spatial = "off",
    pop_spatiotemporal = "off",
    q_diffs_system = "off",
    q_diffs_time = "off",
    q_diffs_spatial = "off",
    index = "monthly",
    extrapolation_grid = grid_common,
    silent = TRUE
  )

  expect_true(fit$settings$extrapolation_grid_supplied)
  expect_equal(fit$data_tmb$n_g, nrow(grid_common))
})
