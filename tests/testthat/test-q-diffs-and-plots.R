test_that("get_q_diffs returns centered temporal effects and projected spatial effects", {
  fit <- mock_jointCPUE_fit()
  qd <- get_q_diffs(
    fit,
    fleet_labels = c("Reference", "Fleet A", "Fleet B"),
    time_values = c(2001, 2002, 2004)
  )

  expect_equal(qd$system$effect, c(0, 0.4, -0.2))
  expect_equal(qd$system$fleet, c("Reference", "Fleet A", "Fleet B"))

  fleet_a <- qd$time[qd$time$fleet == "Fleet A", ]
  fleet_b <- qd$time[qd$time$fleet == "Fleet B", ]

  expect_equal(fleet_a$time, c(2001, 2002, 2004))
  expect_equal(fleet_a$effect, c(-0.5, 0.5, NA))
  expect_equal(fleet_a$observed, c(TRUE, TRUE, FALSE))
  expect_equal(fleet_b$effect, c(-10, NA, 10))
  expect_equal(fleet_b$observed, c(TRUE, FALSE, TRUE))

  expect_equal(qd$spatial$effect, c(0.1, 0.2, -0.1, -0.2))
  expect_equal(unique(qd$spatial$fleet), c("Fleet A", "Fleet B"))
})

test_that("plot helpers return ggplot objects for available q-differences and index input", {
  fit <- mock_jointCPUE_fit()
  index_df <- data.frame(
    time = 1:3,
    index = c(10, 12, 15),
    cv = c(0.2, 0.15, 0.1)
  )

  expect_s3_class(plot_q_diffs_system(fit), "ggplot")
  expect_s3_class(plot_q_diffs_time(fit, time_values = c(2001, 2002, 2004)), "ggplot")
  expect_s3_class(plot_q_diffs_spatial(fit), "ggplot")
  month_plot <- plot_month_diffs(fit)
  expect_s3_class(month_plot, "ggplot")
  expect_s3_class(
    plot_index(index_df, time_values = c("2001", "2002", "2004"), time_positions = c(1, 2, 4)),
    "ggplot"
  )
})

test_that("plot_index auto-rotates year-month labels vertically", {
  index_df <- data.frame(
    time = 1:3,
    index = c(10, 12, 15),
    cv = c(0.2, 0.15, 0.1)
  )

  p_month <- plot_index(index_df, time_values = c("2001-01", "2001-02", "2001-04"))
  p_year <- plot_index(index_df, time_values = c("2001", "2002", "2004"))

  expect_equal(p_month$theme$axis.text.x$angle, 90)
  expect_equal(p_month$theme$axis.text.x$vjust, 0.5)
  expect_equal(p_month$labels$x, "Year-Month")
  expect_equal(p_year$theme$axis.text.x$angle, 45)
  expect_equal(p_year$labels$x, "Year")
})

test_that("get_month_diffs returns observed month values and effects", {
  fit <- mock_jointCPUE_fit(month_values = c(4L, 6L, 11L))
  md <- get_month_diffs(fit)

  expect_equal(md$month, c(4L, 6L, 11L))
  expect_equal(md$effect, c(0.1, 0.2, -0.3))
})

test_that("plot and extraction helpers fail clearly when required components are absent", {
  fit <- mock_jointCPUE_fit(
    q_diffs_system = "off",
    q_diffs_time = "off",
    q_diffs_spatial = "off",
    month_diffs = "off"
  )

  expect_null(get_q_diffs(fit)$system)
  expect_null(get_q_diffs(fit)$time)
  expect_null(get_q_diffs(fit)$spatial)
  expect_null(get_month_diffs(fit))
  expect_error(plot_q_diffs_system(fit), "Systematic q-differences are not available")
  expect_error(plot_q_diffs_time(fit), "Temporal q-differences are not available")
  expect_error(plot_q_diffs_spatial(fit), "Spatial q-differences are not available")
  expect_error(plot_month_diffs(fit), "Month fixed effects are not available")
  expect_error(plot_index(data.frame(time = 1:3, index = 1:3)), "`index_df` is missing required columns")
})
