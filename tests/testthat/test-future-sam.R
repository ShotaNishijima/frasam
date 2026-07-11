context("future SAM workflow")

test_that("simulated SAM results can be combined and projected", {
  result <- make_short_simulated_future_data(n_sim = 2)

  expect_equal(dim(result$par_sim)[2], 2)
  expect_true(all(vapply(result$sam_sim, inherits, logical(1), "sam")))
  expect_true(all(vapply(
    result$sam_sim,
    function(x) x$opt$convergence == 0,
    logical(1)
  )))

  expect_equal(length(result$future_list), 2)
  expect_equal(result$future$data$nsim, 2)
  expect_equal(dim(result$future$data$faa_mat)[3], 2)

  faa_last <- result$future$data$faa_mat[
    , dim(result$future$data$faa_mat)[2], , drop = FALSE
  ]
  expect_true(all(is.finite(faa_last)))
  expect_true(all(faa_last > 0))

  projection <- frasyr::future_vpa(
    tmb_data = result$future$data,
    optim_method = "none",
    multi_init = 1
  )

  expect_s3_class(projection, "future_new")
  expect_equal(dim(projection$faa)[3], 2)
  expect_true(all(is.finite(projection$faa)))
})