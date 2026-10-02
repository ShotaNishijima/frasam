test_that("VPA without sdreport can be plotted without confidence intervals", {
  load(testthat::test_path("..", "..", "data", "res_vpa.rda"))
  res_vpa$rep <- NULL

  plot <- plot_samvpa(res_vpa, CI = 0, scenario_name = "VPA",
                      years = 2010:2013)
  expect_s3_class(plot, "ggplot")
  expect_true(all(plot$data$Year %in% 2010:2013))
  expect_equal(levels(plot$data$Model), "VPA")
  expect_false(any(c("CV", "lower", "upper") %in% names(plot$data)))
  expect_no_error(ggplot2::ggplot_build(plot))
  expect_error(plot_samvpa(res_vpa, CI = 0.95),
               "Rerun vpa\\(\\) with TMB=TRUE & sdreport=TRUE!")
})

test_that("SAM confidence intervals remain available", {
  load(testthat::test_path("..", "..", "data", "res_bh.rda"))
  plot <- plot_samvpa(res_bh, CI = 0.95)
  expect_true(all(c("CV", "lower", "upper") %in% names(plot$data)))
  expect_no_error(ggplot2::ggplot_build(plot))

  res_bh$rep <- NULL
  expect_no_error(ggplot2::ggplot_build(plot_samvpa(res_bh, CI = 0)))
})
