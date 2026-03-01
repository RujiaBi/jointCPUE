test_that("make_data rejects non-positive cpue values", {
  data_input <- tiny_cpue_data(cpue = c(0, 2, 3, 4))
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )

  expect_error(
    make_data(utm$data_utm, mesh),
    "`cpue` must be strictly positive and finite"
  )
})

test_that("make_data yearly mode recodes observed month levels without requiring all 12 months", {
  data_input <- tiny_cpue_data(month = c(11, 4, 11, 6))
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )

  prep <- make_data(utm$data_utm, mesh, index = "yearly")

  expect_equal(prep$data$n_m, 3L)
  expect_equal(prep$data$month_i, c(2L, 0L, 2L, 1L))
  expect_equal(prep$data$use_month_fe, 1L)
  expect_equal(prep$time$month_diffs, "on")
  expect_equal(prep$time$month_values, c(4L, 6L, 11L))
  expect_equal(prep$time$month_col, "month")
})

test_that("make_data yearly mode requires a valid calendar month column", {
  data_input <- tiny_cpue_data()
  data_input$month <- NULL
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )

  expect_error(
    make_data(utm$data_utm, mesh, index = "yearly"),
    "requires a calendar month column"
  )
})

test_that("make_data yearly mode can disable month fixed effects", {
  data_input <- tiny_cpue_data()
  data_input$month <- NULL
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )

  prep <- make_data(
    utm$data_utm,
    mesh,
    index = "yearly",
    month_diffs = "off"
  )

  expect_equal(prep$data$n_m, 1L)
  expect_equal(prep$data$month_i, rep.int(0L, nrow(utm$data_utm)))
  expect_equal(prep$data$use_month_fe, 0L)
  expect_equal(prep$time$month_diffs, "off")
  expect_null(prep$time$month_col)
  expect_equal(prep$time$month_values, integer(0))
})

test_that("make_data monthly mode keeps month fixed effects off", {
  data_input <- tiny_cpue_data(month = c(1, 2, 3, 4))
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )

  prep <- make_data(
    utm$data_utm,
    mesh,
    index = "monthly",
    month_diffs = "on"
  )

  expect_equal(prep$data$n_m, 1L)
  expect_equal(prep$data$use_month_fe, 0L)
  expect_equal(prep$time$month_diffs, "off")
  expect_null(prep$time$month_col)
})

test_that("make_data can reuse a supplied raw extrapolation grid", {
  data_input <- tiny_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )

  grid_common <- unique(utm$data_utm[, c("lon", "lat", "utm_x_scale", "utm_y_scale")])

  prep_auto <- make_data(utm$data_utm, mesh)
  prep_grid <- make_data(utm$data_utm, mesh, extrapolation_grid = grid_common)

  expect_equal(prep_grid$data$n_g, nrow(grid_common))
  expect_equal(prep_grid$key$area_km2, prep_auto$key$area_km2)
  expect_equal(dim(prep_grid$data$A_gs), dim(prep_auto$data$A_gs))
})

test_that("make_data can reuse a supplied extrapolation key with precomputed area", {
  data_input <- tiny_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )

  prep_auto <- make_data(utm$data_utm, mesh)

  extra_key <- data.frame(
    utm_x_scale = mean(range(prep_auto$key$utm_x_scale)) + 0.01,
    utm_y_scale = mean(range(prep_auto$key$utm_y_scale)) + 0.01,
    area_km2 = prep_auto$key$area_km2[1]
  )

  grid_common <- rbind(
    prep_auto$key[, c("utm_x_scale", "utm_y_scale", "area_km2")],
    extra_key
  )

  prep_grid <- make_data(utm$data_utm, mesh, extrapolation_grid = grid_common)

  expect_equal(prep_grid$data$n_g, nrow(grid_common))
  expect_equal(length(prep_grid$data$area_g), nrow(grid_common))
  expect_equal(dim(prep_grid$data$A_gs), c(nrow(grid_common), ncol(prep_auto$data$A_gs)))
})
