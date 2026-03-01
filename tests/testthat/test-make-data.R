tiny_cpue_data <- function(cpue = c(1, 2, 3, 4)) {
  data.frame(
    cpue = cpue,
    lon = c(150, 150.5, 151, 151.5),
    lat = c(40, 40.5, 41, 41.5),
    tid = c(0, 0, 1, 1),
    fleetid = c(0, 1, 0, 1)
  )
}

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
