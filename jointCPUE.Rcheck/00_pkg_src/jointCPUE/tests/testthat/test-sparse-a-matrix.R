test_that("Ais_ij and Ais_x come from the same Tsparse representation", {
  data_input <- tiny_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )
  prep <- make_data(utm$data_utm, mesh)

  A_isT <- methods::as(prep$data$A_is, "TsparseMatrix")

  expect_equal(prep$data$Ais_ij, cbind(A_isT@i, A_isT@j))
  expect_equal(prep$data$Ais_x, A_isT@x)
})
