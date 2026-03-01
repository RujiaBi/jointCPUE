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

test_that("supplying the same extrapolation grid reproduces auto-grid results", {
  data_input <- data.frame(
    cpue = c(1.1, 2.0, 1.7, 2.8, 1.3, 2.4),
    lon = c(150.0, 150.5, 151.0, 150.1, 150.6, 151.1),
    lat = c(40.0, 40.5, 41.0, 40.1, 40.6, 41.1),
    tid = c(0, 0, 0, 1, 1, 1),
    fleetid = c(0, 1, 0, 0, 1, 0)
  )
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )
  grid_common <- unique(utm$data_utm[, c("lon", "lat", "utm_x_scale", "utm_y_scale")])

  prep_auto <- make_data(utm$data_utm, mesh, index = "monthly")
  prep_grid <- make_data(utm$data_utm, mesh, index = "monthly", extrapolation_grid = grid_common)

  expect_equal(prep_auto$key$area_km2, prep_grid$key$area_km2)
  expect_equal(prep_auto$key$area_km2_scaled, prep_grid$key$area_km2_scaled)
  expect_equal(as.matrix(prep_auto$data$A_gs), as.matrix(prep_grid$data$A_gs))

  fit_auto <- jointCPUE(
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

  fit_grid <- jointCPUE(
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

  index_auto <- get_index(fit_auto)
  index_grid <- get_index(fit_grid)

  expect_equal(index_auto$index, index_grid$index)
  expect_equal(index_auto$cv, index_grid$cv)
})
