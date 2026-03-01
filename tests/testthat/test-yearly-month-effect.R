test_that("yearly fits can include an observed-month fixed effect without changing index length", {
  data_input <- data.frame(
    cpue = c(1.0, 1.2, 1.4, 2.0, 2.2, 2.4),
    lon = c(150.0, 150.2, 150.4, 150.1, 150.3, 150.5),
    lat = c(40.0, 40.1, 40.2, 40.3, 40.4, 40.5),
    tid = c(0, 0, 0, 1, 1, 1),
    fleetid = rep(0, 6),
    month = c(4, 6, 11, 4, 6, 11)
  )

  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.1
  )

  fit <- jointCPUE(
    data_utm = utm$data_utm,
    mesh = mesh,
    pop_spatial = "off",
    pop_spatiotemporal = "off",
    q_diffs_system = "off",
    q_diffs_time = "off",
    q_diffs_spatial = "off",
    index = "yearly",
    month_diffs = "on",
    silent = TRUE
  )

  idx <- get_index(fit)

  expect_equal(fit$settings$index, "yearly")
  expect_equal(fit$settings$month_diffs, "on")
  expect_equal(fit$settings$month_values, c(4L, 6L, 11L))
  expect_equal(fit$data_tmb$use_month_fe, 1L)
  expect_equal(nrow(idx), 2L)
  expect_true(all(is.finite(idx$index)))
})

test_that("yearly fits can disable month fixed effects", {
  data_input <- data.frame(
    cpue = c(1.0, 1.2, 1.4, 2.0, 2.2, 2.4),
    lon = c(150.0, 150.2, 150.4, 150.1, 150.3, 150.5),
    lat = c(40.0, 40.1, 40.2, 40.3, 40.4, 40.5),
    tid = c(0, 0, 0, 1, 1, 1),
    fleetid = rep(0, 6)
  )

  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.1
  )

  fit <- jointCPUE(
    data_utm = utm$data_utm,
    mesh = mesh,
    pop_spatial = "off",
    pop_spatiotemporal = "off",
    q_diffs_system = "off",
    q_diffs_time = "off",
    q_diffs_spatial = "off",
    index = "yearly",
    month_diffs = "off",
    silent = TRUE
  )

  idx <- get_index(fit)

  expect_equal(fit$settings$index, "yearly")
  expect_equal(fit$settings$month_diffs, "off")
  expect_equal(fit$settings$month_values, integer(0))
  expect_equal(fit$data_tmb$use_month_fe, 0L)
  expect_equal(nrow(idx), 2L)
  expect_true(all(is.finite(idx$index)))
})

test_that("monthly fits always keep month fixed effects off", {
  data_input <- data.frame(
    cpue = c(1.0, 1.2, 1.4, 2.0),
    lon = c(150.0, 150.2, 150.4, 150.1),
    lat = c(40.0, 40.1, 40.2, 40.3),
    tid = c(0, 1, 2, 3),
    fleetid = rep(0, 4),
    month = c(4, 6, 8, 11)
  )

  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.1
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
    month_diffs = "on",
    silent = TRUE
  )

  expect_equal(fit$settings$index, "monthly")
  expect_equal(fit$settings$month_diffs, "off")
  expect_equal(fit$data_tmb$use_month_fe, 0L)
})
