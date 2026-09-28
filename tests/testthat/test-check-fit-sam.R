library(frasam)

context("check_fit_sam")

test_that("check_fit_sam returns convergence diagnostics", {
  data("samres_example", package = "frasam")

  fit_check <- check_fit_sam(samres, verbose = FALSE)

  expect_true(is.list(fit_check))
  expect_true(is.logical(fit_check$ok))
  expect_true(is.data.frame(fit_check$checks))
  expect_true(all(c("check", "ok", "value", "threshold", "message") %in% names(fit_check$checks)))
  expect_true(is.data.frame(fit_check$fixed))
  expect_true(is.data.frame(fit_check$sigma))
  expect_true(all(c("type", "index", "value", "ok", "problem") %in% names(fit_check$sigma)))
})

test_that("check_fit_sam detects a failed gradient diagnostic", {
  data("samres_example", package = "frasam")

  fit_check <- check_fit_sam(samres, gradient_tol = 0, verbose = FALSE)

  gradient_row <- fit_check$checks[fit_check$checks$check == "maximum fixed-effect gradient", ]
  expect_false(gradient_row$ok)
})

test_that("check_fit_sam reports which sigma failed the range check", {
  data("samres_example", package = "frasam")

  fit_check <- check_fit_sam(samres, sigma_range = c(0.5, 1), verbose = FALSE)

  sigma_row <- fit_check$checks[fit_check$checks$check == "reported sigma range", ]
  expect_false(sigma_row$ok)
  expect_true(any(!fit_check$sigma$ok))
  expect_true(all(c("too small", "too large") %in% unique(fit_check$sigma$problem[!fit_check$sigma$ok])))
})
