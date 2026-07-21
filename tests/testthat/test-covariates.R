make_covariate_cpue_data <- function() {
  data.frame(
    cpue = c(1.1, 2.0, 1.7, 2.8, 1.3, 2.4, 1.9, 2.2, 1.5, 2.6, 1.8, 2.1),
    lon = c(150.0, 150.1, 150.2, 150.3, 150.4, 150.5, 150.6, 150.7, 150.8, 150.9, 151.0, 151.1),
    lat = c(40.0, 40.0, 40.1, 40.1, 40.2, 40.2, 40.3, 40.3, 40.4, 40.4, 40.5, 40.5),
    tid = rep(0:2, each = 4),
    fleetid = rep(c(0L, 1L), 6),
    month = rep(4:9, length.out = 12),
    depth = c(100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210),
    sst = seq(10, 21, length.out = 12),
    season = rep(c("spring", "summer", "fall"), each = 4)
  )
}

test_that("make_data parses catchability and population smooths", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("mgcv")

  data_input <- make_covariate_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  prep <- make_data(
    data_utm = utm$data_utm,
    mesh = mesh,
    formula_catchability = ~ s(depth),
    formula_population = ~ s(sst),
    index = "monthly"
  )

  expect_equal(prep$data$has_smooths_catch, 1L)
  expect_equal(prep$data$has_smooths_pop, 1L)
  expect_equal(nrow(prep$data$Xs_catch), nrow(utm$data_utm))
  expect_equal(nrow(prep$data$Xs_pop_i), nrow(utm$data_utm))
  expect_equal(nrow(prep$data$Xs_pop_g), prep$data$n_g * prep$data$n_t)
  expect_equal(length(prep$data$Zs_pop_i), length(prep$data$Zs_pop_g))
  expect_equal(nrow(prep$projection_data), prep$data$n_g * prep$data$n_t)
})

test_that("fixed catchability factors use sum-to-zero coding", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")

  data_input <- make_covariate_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  prep <- make_data(
    data_utm = utm$data_utm,
    mesh = mesh,
    formula_catchability = ~ season,
    index = "monthly"
  )

  basis <- prep$covariate_info$catchability$fixed_terms[[1]]
  beta <- c(0.2, -0.1)
  effects <- as.numeric(prep$data$Xf_catch %*% beta)
  means <- tapply(effects, utm$data_utm$season, unique)

  expect_equal(basis$kind, "factor")
  expect_equal(ncol(prep$data$Xf_catch), length(basis$levels) - 1L)
  expect_equal(unname(sum(unlist(means))), 0, tolerance = 1e-12)
})

test_that("random factor covariates use f(var) syntax", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")

  data_input <- make_covariate_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  prep <- make_data(
    data_utm = utm$data_utm,
    mesh = mesh,
    formula_catchability = ~ f(season),
    formula_population = ~ f(season),
    index = "monthly"
  )

  catch_basis <- prep$covariate_info$catchability$random_factor_terms[[1]]
  pop_basis <- prep$covariate_info$population$random_factor_terms[[1]]

  expect_equal(prep$data$has_random_factors_catch, 1L)
  expect_equal(prep$data$has_random_factors_pop, 1L)
  expect_equal(catch_basis$kind, "random_factor")
  expect_equal(pop_basis$kind, "random_factor")
  expect_equal(length(prep$data$Zf_catch), 1L)
  expect_equal(length(prep$data$Zf_pop_i), 1L)
  expect_equal(length(prep$data$Zf_pop_g), 1L)
  expect_equal(nrow(prep$data$Zf_catch[[1]]), nrow(utm$data_utm))
  expect_equal(nrow(prep$data$Zf_pop_g[[1]]), prep$data$n_g * prep$data$n_t)
  expect_equal(ncol(prep$data$Zf_catch[[1]]), length(catch_basis$levels))
})

test_that("population covariates require projection data when varying within a cell", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("mgcv")

  data_input <- make_covariate_cpue_data()
  data_input$lon[2] <- data_input$lon[1]
  data_input$lat[2] <- data_input$lat[1]
  data_input$sst[2] <- data_input$sst[1] + 5
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  expect_error(
    make_data(
      data_utm = utm$data_utm,
      mesh = mesh,
      formula_population = ~ s(sst),
      index = "monthly"
    ),
    "varies within an extrapolation grid cell"
  )

  prep0 <- make_data(utm$data_utm, mesh, index = "monthly")
  projection_data <- prep0$key[, c("utm_x_scale", "utm_y_scale"), drop = FALSE]
  projection_data$sst <- seq(30, 30 + nrow(projection_data) - 1)

  prep <- make_data(
    data_utm = utm$data_utm,
    mesh = mesh,
    formula_population = ~ s(sst),
    projection_data = projection_data,
    index = "monthly"
  )

  expect_equal(nrow(prep$data$Xs_pop_g), prep$data$n_g * prep$data$n_t)
})

test_that("month_diffs prevents duplicate month covariates in formulas", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")

  data_input <- make_covariate_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.2
  )

  expect_error(
    make_data(
      data_utm = utm$data_utm,
      mesh = mesh,
      index = "yearly",
      month_diffs = "on",
      formula_catchability = ~ month
    ),
    "already adds `month`"
  )
})
