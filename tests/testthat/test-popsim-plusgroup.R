test_that("SAM population simulations use aggregate catch errors and preserve NA", {
  predicted <- matrix(c(10, 20, 30, 12, 24, 36), 3,
                      dimnames = list(0:2, 2000:2001))
  observed <- predicted
  observed[2, 1] <- 50
  observed[3, 1] <- NA_real_
  obs <- cbind(year = rep(2000:2001, each = 3), fleet = 1,
               age = rep(0:2, 2), obs = as.vector(observed),
               maxage = rep(0:2, 2))
  res <- structure(list(
    input = list(dat = list(caa = as.data.frame(observed),
                            index = matrix(c(1, NA_real_), 1)),
                 change_plusgroup = TRUE, last.catch.zero = FALSE),
    caa = predicted, pred.index = matrix(1, 1, 2),
    sigma = 0.1, sigma.logC = c(0.2, 0.3, 0.5),
    data = list(obs = fix_maxage_and_remove_na(obs, verbose = FALSE))
  ), class = "sam")
  simulate <- popsim_vpasam
  environment(simulate) <- new.env(parent = environment(popsim_vpasam))
  calls <- list()
  environment(simulate)$rnorm <- function(n, mean, sd) {
    calls[[length(calls) + 1L]] <<- list(n = n, mean = mean, sd = sd)
    mean + sd
  }
  aggregate.sd <- sqrt((20 * 0.3)^2 + (30 * 0.5)^2) / 50
  results <- simulate(res, n = 2, seed = NULL)
  expect_equal(length(results), 2L)
  expect_equal(length(calls), 8L)
  expect_equal(as.numeric(calls[[3]]$mean), log(c(50, 24)))
  expect_equal(as.numeric(calls[[3]]$sd), c(aggregate.sd, 0.3))
  expect_equal(calls[[4]]$n, 1L)
  for (dat in results) {
    expect_equal(dat$caa[2, 1], exp(log(50) + aggregate.sd))
    expect_equal(dat$caa[2, 2], exp(log(24) + 0.3))
    expect_identical(is.na(dat$caa), is.na(res$input$dat$caa))
    expect_identical(is.na(dat$index), is.na(res$input$dat$index))
  }
  # The generated data keep the missing-age marker needed on refitting.
  generated.obs <- obs
  generated.obs[, "obs"] <- as.vector(as.matrix(results[[1]]$caa))
  fixed <- fix_maxage_and_remove_na(generated.obs, verbose = FALSE)
  expect_true(is.matrix(fixed))
  expect_equal(nrow(fixed), 5L)
  expect_equal(unname(fixed[2, "maxage"]), 2)

  res$input$last.catch.zero <- TRUE
  result <- simulate(res, n = 1, seed = NULL)[[1]]
  expect_equal(as.numeric(result$caa[, 2]), rep(0, 3))
  expect_true(is.na(result$caa[3, 1]))

  res$input$last.catch.zero <- FALSE
  res$input$change_plusgroup <- FALSE
  calls <- list()
  result <- simulate(res, n = 1, seed = NULL)[[1]]
  expect_equal(as.numeric(calls[[3]]$mean), log(c(20, 24)))
  expect_equal(as.numeric(calls[[3]]$sd), c(0.3, 0.3))
  expect_equal(result$caa[2, 1], exp(log(20) + 0.3))

  # Real random draws remain reproducible with the supplied seed.
  expect_identical(popsim_vpasam(res, n = 2, seed = 123),
                   popsim_vpasam(res, n = 2, seed = 123))
})
