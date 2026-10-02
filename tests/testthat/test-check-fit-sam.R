library(frasam)

# A controlled fit isolates optimizer status from all other diagnostics.
diagnostic_fit <- function() {
  structure(list(
    opt = list(convergence = 1L, message = "false convergence (8)", par = c(x = 1)),
    rep = list(pdHess = TRUE, gradient.fixed = c(x = 0),
               par.fixed = c(x = 1), cov.fixed = matrix(1, 1, 1)),
    input = list()
  ), class = "sam")
}

test_that("optimizer status is advisory by default and required in strict mode", {
  fit <- diagnostic_fit()
  result <- check_fit_sam(fit, verbose = FALSE)
  row <- result$checks[result$checks$check == "optimizer convergence", ]
  expect_true(result$ok)
  expect_false(row$ok)
  expect_false(row$required)
  expect_equal(row$message, "false convergence (8)")
  expect_equal(result$ok, all(result$checks$ok[result$checks$required]))
  expect_message(check_fit_sam(fit), "inspect the stopping message")

  strict <- check_fit_sam(fit, require_optimizer_convergence = TRUE, verbose = FALSE)
  expect_false(strict$ok)
  expect_true(all(strict$checks$required))
  fit$opt$convergence <- 0L
  expect_true(check_fit_sam(fit, require_optimizer_convergence = TRUE, verbose = FALSE)$ok)
})

test_that("gradient and Hessian failures still fail the overall result", {
  fit <- diagnostic_fit()
  fit$rep$gradient.fixed <- c(x = 0.1)
  expect_false(check_fit_sam(fit, verbose = FALSE)$ok)
  fit$rep$gradient.fixed <- c(x = 0)
  fit$rep$pdHess <- FALSE
  expect_false(check_fit_sam(fit, verbose = FALSE)$ok)
})

test_that("strict mode requires a single logical value", {
  for (value in list(NA, NULL, c(TRUE, FALSE), "TRUE", 1)) {
    expect_error(check_fit_sam(diagnostic_fit(), require_optimizer_convergence = value),
                 "must be TRUE or FALSE")
  }
})

context("check_fit_sam")

test_that("check_fit_sam returns convergence diagnostics", {
  data("samres_example", package = "frasam")

  fit_check <- check_fit_sam(samres, verbose = FALSE)

  expect_true(is.list(fit_check))
  expect_true(is.logical(fit_check$ok))
  expect_true(is.data.frame(fit_check$checks))
  expect_true(all(c("check", "ok", "value", "threshold", "message", "required") %in% names(fit_check$checks)))
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
