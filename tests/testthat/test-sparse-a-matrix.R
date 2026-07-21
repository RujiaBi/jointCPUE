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

test_that("mesh basis is recomputed when coordinates change without changing row count", {
  data_input <- tiny_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )

  loc_reordered <- mesh$loc_xy[rev(seq_len(nrow(mesh$loc_xy))), , drop = FALSE]
  mesh_updated <- jointCPUE:::.as_jointCPUEmesh(
    mesh,
    loc_xy = loc_reordered,
    recompute_A = "auto"
  )
  expected <- fmesher::fm_basis(mesh$mesh, loc = loc_reordered)

  expect_equal(as.matrix(mesh_updated$A), as.matrix(expected))
  expect_false(isTRUE(all.equal(as.matrix(mesh_updated$A), as.matrix(mesh$A))))
})

test_that("make_data rejects observation and projection coordinates outside the mesh", {
  data_input <- tiny_cpue_data()
  utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
  mesh <- make_mesh(
    utm$data_utm,
    xy_cols = c("utm_x_scale", "utm_y_scale"),
    type = "cutoff",
    cutoff = 0.5
  )

  data_outside <- utm$data_utm
  data_outside$utm_x_scale[1] <- max(mesh$mesh$loc[, 1]) + 100
  data_outside$utm_y_scale[1] <- max(mesh$mesh$loc[, 2]) + 100
  expect_error(
    make_data(data_outside, mesh),
    "`A_is` contains invalid FEM projection rows"
  )

  prep <- make_data(utm$data_utm, mesh)
  grid_outside <- rbind(
    prep$key[, c("utm_x_scale", "utm_y_scale", "area_km2")],
    data.frame(
      utm_x_scale = max(mesh$mesh$loc[, 1]) + 100,
      utm_y_scale = max(mesh$mesh$loc[, 2]) + 100,
      area_km2 = prep$key$area_km2[1]
    )
  )
  expect_error(
    make_data(utm$data_utm, mesh, extrapolation_grid = grid_outside),
    "`A_gs` contains invalid FEM projection rows"
  )
})
