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
