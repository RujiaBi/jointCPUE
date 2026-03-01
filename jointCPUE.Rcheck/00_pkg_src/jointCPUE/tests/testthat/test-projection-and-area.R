test_that("make_utm normalizes 0..360 longitudes and adds scaled coordinates", {
  data_input <- data.frame(
    lon = c(359, 1),
    lat = c(10, 10)
  )

  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")

  expect_equal(utm$data_utm$lon_std, c(-1, 1))
  expect_equal(utm$utm_zone, 31L)
  expect_true(all(c("utm_x", "utm_y", "utm_x_scale", "utm_y_scale") %in% names(utm$data_utm)))
  expect_true(all(is.finite(utm$data_utm$utm_x_scale)))
  expect_true(all(is.finite(utm$data_utm$utm_y_scale)))
})

test_that("calculate_area returns positive areas aligned to input order", {
  grid <- data.frame(
    X = c(0, 1, 0, 1),
    Y = c(10, 10, 11, 11),
    utm_x_scale = c(1, 2, 1, 2),
    utm_y_scale = c(3, 3, 4, 4)
  )

  area <- calculate_area(grid, cellsize = c(1, 1), check_grid = FALSE)

  expect_equal(nrow(area), nrow(grid))
  expect_equal(area$order, seq_len(nrow(grid)))
  expect_true(all(c("utm_x_scale", "utm_y_scale", "area_km2") %in% names(area)))
  expect_true(all(is.finite(area$area_km2)))
  expect_true(all(area$area_km2 > 0))
})
