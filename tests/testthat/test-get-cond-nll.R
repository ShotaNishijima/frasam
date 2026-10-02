library(frasam)

test_that("get_cond_nll aggregates conditional NLL by component and fleet", {
  report <- list(ans_n = 10, ans_f = 20, ans_obs = c(1, 2, 3, 4, 5))
  samres <- structure(
    list(
      obj = list(report = function() report),
      data = list(obs = data.frame(fleet = c(1, 1, 2, 3, 3)))
    ),
    class = "sam"
  )

  out <- get_cond_nll(samres, index_name = c("Survey A", "Survey B"))

  expect_s3_class(out, "tbl_df")
  expect_identical(
    as.character(out$type),
    c("Process_N", "Process_F", "Catch_at_age", "Survey A", "Survey B")
  )
  expect_equal(out$nll, c(10, 20, 3, 3, 9))
  expect_identical(levels(out$type), as.character(out$type))
})

test_that("get_cond_nll validates the TMB report and index names", {
  samres <- structure(
    list(
      obj = list(report = function() list(ans_n = 1, ans_f = 2)),
      data = list(obs = data.frame(fleet = c(1, 2)))
    ),
    class = "sam"
  )

  expect_error(get_cond_nll(samres), "missing: ans_obs", fixed = TRUE)

  samres$obj$report <- function() {
    list(ans_n = 1, ans_f = 2, ans_obs = c(3, 4))
  }
  expect_error(get_cond_nll(samres, index_name = character()), "index_name")
})
