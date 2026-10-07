test_that("plus-group observations always remain a matrix", {
  obs <- cbind(year = c(2000, 2000), fleet = 1,
               age = 0:1, obs = c(10, 20), maxage = 0:1)
  expect_identical(fix_maxage_and_remove_na(obs, verbose = FALSE), obs)
  obs[2, "obs"] <- NA_real_
  fixed <- fix_maxage_and_remove_na(obs, verbose = FALSE)
  expect_true(is.matrix(fixed))
  expect_equal(unname(fixed[1, "maxage"]), 1)
  expect_equal(nrow(fixed), 1L)
})

test_that("catch bootstraps preserve missing ages and aggregate plus groups", {
  predicted <- matrix(c(10, 20, 30, 10, 20, 30), 3,
                      dimnames = list(0:2, 2000:2001))
  observed <- predicted
  observed[2, 1] <- 50
  observed[3, 1] <- NA_real_
  obs <- cbind(year = rep(2000:2001, each = 3), fleet = 1,
               age = rep(0:2, 2), obs = as.vector(observed),
               maxage = rep(0:2, 2))
  res <- list(
    input = list(dat = list(caa = as.data.frame(observed),
                            index = matrix(1, 1, 2)),
                 change_plusgroup = TRUE, last.catch.zero = FALSE),
    caa = predicted, pred.index = matrix(1, 1, 2),
    sigma = 0, sigma.logC = rep(0, 3),
    data = list(obs = fix_maxage_and_remove_na(obs, verbose = FALSE))
  )
  # Capture the simulated input without fitting a different model in each test.
  boot <- boo_sam
  environment(boot) <- new.env(parent = environment(boo_sam))
  environment(boot)$sam <- function(dat, ...) list(dat = dat)
  for (method in c("p", "n")) {
    result <- boot(res, n = 1, method = method)[[1]]$dat$caa
    expect_identical(is.na(result), is.na(res$input$dat$caa))
    expect_equal(as.matrix(result), observed)
  }
  # A single available residual must be sampled by position, even when nonzero.
  res$input$dat$caa[3, 2] <- 60
  result <- boot(res, n = 1, method = "n")[[1]]$dat$caa
  expect_equal(result[3, 2], 60)
})

test_that("parametric catch bootstrap uses the aggregate mean and sigma in rnorm", {
  predicted <- matrix(c(10, 20, 30, 12, 24, 36), 3,
                      dimnames = list(0:2, 2000:2001))
  observed <- predicted
  observed[2, 1] <- 50
  observed[3, 1] <- NA_real_
  obs <- cbind(year = rep(2000:2001, each = 3), fleet = 1,
               age = rep(0:2, 2), obs = as.vector(observed),
               maxage = rep(0:2, 2))
  res <- list(
    input = list(dat = list(caa = as.data.frame(observed),
                            index = matrix(1, 1, 2)),
                 change_plusgroup = TRUE, last.catch.zero = FALSE),
    caa = predicted, pred.index = matrix(1, 1, 2),
    sigma = 0.1, sigma.logC = c(0.2, 0.3, 0.5),
    data = list(obs = fix_maxage_and_remove_na(obs, verbose = FALSE))
  )
  boot <- boo_sam
  environment(boot) <- new.env(parent = environment(boo_sam))
  environment(boot)$sam <- function(dat, ...) list(dat = dat)
  calls <- list()
  # A fixed standard-normal draw makes the generated values deterministic,
  # while capturing the actual mean and sd supplied by boo_sam.
  environment(boot)$rnorm <- function(n, mean, sd) {
    calls[[length(calls) + 1L]] <<- list(n = n, mean = mean, sd = sd)
    mean + sd
  }
  result <- boot(res, n = 1, method = "p")[[1]]$dat$caa
  aggregate.sd <- sqrt((20 * 0.3)^2 + (30 * 0.5)^2) / 50

  expect_equal(length(calls), 4L) # One index and three catch-age calls.
  expect_equal(calls[[3]]$n, 2L)
  expect_equal(as.numeric(calls[[3]]$mean), log(c(50, 24)))
  expect_equal(as.numeric(calls[[3]]$sd), c(aggregate.sd, 0.3))
  expect_equal(result[2, 1], exp(log(50) + aggregate.sd))
  expect_equal(result[2, 2], exp(log(24) + 0.3))
  expect_equal(result[1, 1], exp(log(10) + 0.2))
  expect_equal(calls[[4]]$n, 1L) # No draw for the missing oldest age.
  expect_equal(as.numeric(calls[[4]]$mean), log(36))
  expect_equal(as.numeric(calls[[4]]$sd), 0.5)
  expect_identical(is.na(result), is.na(res$input$dat$caa))

  # With the option disabled, use the original age-specific mean and sigma.
  res$input$change_plusgroup <- FALSE
  calls <- list()
  result <- boot(res, n = 1, method = "p")[[1]]$dat$caa
  expect_equal(as.numeric(calls[[3]]$mean), log(c(20, 24)))
  expect_equal(as.numeric(calls[[3]]$sd), c(0.3, 0.3))
  expect_equal(result[2, 1], exp(log(20) + 0.3))
})
