library(testthat)
library(dplyr)
library(tibble)

# Isolate candidate fitting while using the real convergence diagnostics.
sigma_grid_env <- new.env(parent = globalenv())
sys.source(testthat::test_path("..", "..", "R", "model_selection.R"), sigma_grid_env)
sys.source(testthat::test_path("..", "..", "R", "check_fit_sam.R"), sigma_grid_env)
grid_fit <- function(aic, gradient = 0, convergence = 0L) {
  structure(list(aic = aic, input = list(),
    opt = list(convergence = convergence, par = c(x = 1)),
    rep = list(pdHess = TRUE, gradient.fixed = c(x = gradient),
               par.fixed = c(x = 1), cov.fixed = matrix(1, 1, 1))), class = "sam")
}

test_that("each stage continues from the best passing candidate", {
  parents <- list()
  sigma_grid_env$divide_sigma <- function(samres, var, which) {
    parents[[length(parents) + 1L]] <<- samres$aic
    if (samres$aic == 100) {
      if (which == 1) grid_fit(10, gradient = 1) else grid_fit(90, convergence = 1L)
    } else grid_fit(80)
  }
  result <- suppressMessages(sigma_grid_env$select_sigma_grid(
    grid_fit(100), data.frame(var = "varC", X = 1:2)))
  expect_equal(unlist(parents), c(100, 100, 90))
  expect_equal(result$bestres$aic, 80)
  expect_identical(result$tbl_sigma$fit_ok, c(TRUE, FALSE, TRUE, TRUE))
  expect_true(all(result$tbl_sigma$fit_ok[!is.na(result$tbl_sigma$model)]))
})

test_that("no passing candidate stops even when stopAIC is FALSE", {
  sigma_grid_env$divide_sigma <- function(...) grid_fit(10, gradient = 1)
  expect_warning(result <- suppressMessages(sigma_grid_env$select_sigma_grid(
    grid_fit(100), data.frame(var = "varC", X = 1:2), stopAIC = FALSE)), "No eligible models")
  expect_length(result$reslist, 1)
  expect_equal(result$bestres$aic, 100)
  expect_true(all(is.na(result$tbl_sigma$model[-1])))
  expect_warning(result <- suppressMessages(sigma_grid_env$select_sigma_grid(
    grid_fit(1, gradient = 1), data.frame(var = "varC", X = 1))), "No eligible models")
  expect_null(result$bestres)
  expect_true(all(is.na(result$tbl_sigma$model)))
})

test_that("invalid initial fits, ties, nonfinite AIC and disabled checks work", {
  sigma_grid_env$divide_sigma <- function(samres, var, which) grid_fit(20)
  result <- suppressMessages(sigma_grid_env$select_sigma_grid(
    grid_fit(1, gradient = 1), data.frame(var = "varC", X = 1:2)))
  expect_equal(result$bestres$aic, 20)
  expect_equal(sum(result$tbl_sigma$model == "best", na.rm = TRUE), 1)
  expect_length(result$reslist, 2)

  sigma_grid_env$divide_sigma <- function(samres, var, which) {
    if (which == 1) grid_fit(NA_real_) else grid_fit(10, gradient = 1)
  }
  expect_warning(result <- suppressMessages(sigma_grid_env$select_sigma_grid(
    grid_fit(100), data.frame(var = "varC", X = 1:2),
    check_converge = FALSE, stopAIC = TRUE)), "No eligible models")
  expect_equal(result$bestres$aic, 10)

  result <- suppressMessages(sigma_grid_env$select_sigma_grid(
    grid_fit(100), data.frame(var = character(), X = numeric())))
  expect_equal(result$bestres$aic, 100)
  expect_length(result$reslist, 0)
})
