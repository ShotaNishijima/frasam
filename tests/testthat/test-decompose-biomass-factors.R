library(frasam)

context("biomass factor decomposition")

test_that("decompose_biomass_factors returns summaries for SAM biomass", {
  data("samres_example", package = "frasam")

  out <- decompose_biomass_factors(samres, target = "biomass")

  expect_true(is.list(out))
  expect_true(all(c(
    "recruitment", "growth", "fishing", "natural", "process", "terminal_loss",
    "target_quantity", "annual_effect_sum", "annual_effect_sum_with_terminal",
    "annual_change", "residual", "residual_with_terminal", "target",
    "age_aggregated", "percent_by_age", "percent_aggregated", "previous_total"
  ) %in% names(out)))
  expect_identical(
    rownames(out$age_aggregated),
    c("recruitment", "growth", "fishing", "process", "natural", "maturity")
  )
  expect_identical(rownames(out$percent_aggregated), rownames(out$age_aggregated))
  expect_identical(colnames(out$age_aggregated), colnames(out$percent_aggregated))
  expect_length(out$previous_total, ncol(out$age_aggregated))
  expect_true(all(is.na(out$percent_aggregated[, 1])))
})

test_that("decompose_biomass_factors returns summaries for VPA SSB", {
  data("res_vpa_example", package = "frasyr")

  out <- decompose_biomass_factors(res_vpa_example, target = "ssb")

  expect_true(is.list(out))
  expect_identical(
    rownames(out$age_aggregated),
    c("recruitment", "growth", "fishing", "process", "natural", "maturity")
  )
  expect_identical(rownames(out$percent_aggregated), rownames(out$age_aggregated))
  expect_identical(colnames(out$age_aggregated), colnames(out$percent_aggregated))
  expect_true("maturity" %in% names(out$percent_by_age))
  expect_true(is.matrix(out$percent_by_age$maturity))
})

test_that("plot_biomass_factors preserves row order and numeric years", {
  data("samres_example", package = "frasam")

  out <- decompose_biomass_factors(samres, target = "biomass")
  gg <- plot_biomass_factors(out)

  expect_s3_class(gg, "ggplot")
  expect_s3_class(ggplot2::ggplot_build(gg), "ggplot_built")
  expect_identical(levels(gg$data$effect), rownames(out$percent_aggregated))
  expect_true(is.numeric(gg$data$year))
  expect_false(is.factor(gg$data$year))
})

test_that("plot_biomass_factors scales absolute values", {
  data("samres_example", package = "frasam")

  out <- decompose_biomass_factors(samres, target = "biomass")
  unscaled <- plot_biomass_factors(out, type = "absolute", scale = 1)
  scaled <- plot_biomass_factors(out, type = "absolute", scale = 1000)

  expect_equal(scaled$data$value, unscaled$data$value / 1000)
})
