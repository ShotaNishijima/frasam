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

test_that("decompose_biomass_effects harmonizes changed plus groups", {
  naa <- matrix(
    c(
      100, 110, 120,
      50, 55, 60,
      20, 25, 30,
      10, 12, NA
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(0:3, 2000:2002)
  )
  waa <- matrix(
    c(
      1, 1, 1,
      2, 2, 2,
      3, 4, 5,
      6, 8, NA
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = dimnames(naa)
  )
  maa <- matrix(
    c(
      0, 0, 0,
      0.5, 0.5, 0.5,
      0.8, 0.9, 1,
      1, 1, NA
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = dimnames(naa)
  )
  faa <- matrix(
    c(
      0.1, 0.1, 0.1,
      0.2, 0.2, 0.2,
      0.3, 0.4, 0.5,
      0.7, 0.8, NA
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = dimnames(naa)
  )
  M <- matrix(
    c(
      0.4, 0.4, 0.4,
      0.4, 0.4, 0.4,
      0.5, 0.6, 0.7,
      0.8, 0.9, NA
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = dimnames(naa)
  )

  out <- suppressWarnings(decompose_biomass_effects(
    naa = naa, waa = waa, faa = faa, M = M, maa = maa, plus_group = TRUE
  ))

  expected_baa <- colSums(naa[3:4, ] * waa[3:4, ], na.rm = TRUE)
  expected_ssb <- colSums(naa[3:4, ] * waa[3:4, ] * maa[3:4, ], na.rm = TRUE)
  expected_faa_2000 <- -log(sum(naa[3:4, 1] * exp(-faa[3:4, 1])) / sum(naa[3:4, 1]))
  expected_M_2000 <- -log(sum(naa[3:4, 1] * exp(-M[3:4, 1])) / sum(naa[3:4, 1]))

  expect_equal(nrow(out$target_quantity), 3)
  expect_equal(rownames(out$target_quantity), c("0", "1", "2"))
  expect_true(out$plus_group_adjustment$changed)
  expect_equal(out$target_quantity[3, ], expected_ssb)

  collapsed <- suppressWarnings(.harmonize_changed_plus_group(
    naa = naa, waa = waa, faa = faa, M = M, maa = maa
  ))
  expect_equal(collapsed$naa[3, ], colSums(naa[3:4, ], na.rm = TRUE))
  expect_equal(collapsed$naa[3, ] * collapsed$waa[3, ], expected_baa)
  expect_equal(collapsed$naa[3, ] * collapsed$waa[3, ] * collapsed$maa[3, ], expected_ssb)
  expect_equal(collapsed$faa[3, 1], expected_faa_2000)
  expect_equal(collapsed$M[3, 1], expected_M_2000)
})
