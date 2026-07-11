context("slow future SAM workflow")

test_that("est_MSYRP runs for combined simulated SAM results", {
  skip_if_not_slow_tests()
  data_future <- make_short_simulated_future_data(n_sim = 2)$future

  result <- frasyr::est_MSYRP(
    data_future = data_future,
    optim_method = "R",
    candidate_PGY = -1,
    candidate_B0 = -1,
    candidate_Babs = -1,
    calc_yieldcurve = FALSE,
    trace.multi = c(0, 0.25, 0.5, 1, 2, 5, 10)
  )

  expect_true(is.list(result))
  expect_true(all(c("summary", "trace", "res_future_MSY") %in% names(result)))
  expect_true(is.finite(result$res_future_MSY$multi))
  expect_gt(result$res_future_MSY$multi, 0)
})