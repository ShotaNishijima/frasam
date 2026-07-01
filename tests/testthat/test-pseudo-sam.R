library(frasam)

context("pseudo SAM example")

test_that("pseudo SAM example runs and matches saved result", {
  data("sam_ex", package = "frasam")

  expect_error(res_test <- fit_pseudo_sam_example(), NA)

  expect_true(inherits(res_test, "sam"))
  expect_equal(res_test$opt$convergence, sam_ex$opt$convergence)

  testcontents <- c(
    "loglik", "aic", "q", "b",
    "sigma", "sigma.logC", "sigma.logFsta", "rho",
    "F", "faa", "N", "naa", "baa", "ssb", "pred.index"
  )

  for (i in seq_along(testcontents)) {
    expect_equal(
      eval(parse(text = paste0("sam_ex$", testcontents[i]))),
      eval(parse(text = paste0("res_test$", testcontents[i]))),
      tolerance = 1e-3
    )
  }

  expect_equal(sam_ex$opt$par, res_test$opt$par, tolerance = 5e-2)

  expect_error(
    plot_res <- plot_samvpa(
      res_test,
      CI = 0,
      scenario_name = "pseudo"
    ),
    NA
  )
  expect_false(is.null(plot_res))

  expect_error(
    plot_res_scaled <- plot_samvpa(
      res_test,
      CI = 0,
      scenario_name = "pseudo",
      what.plot = c("biomass", "SSB", "Recruitment"),
      scale_biomass = 1,
      scale_ssb = 1,
      scale_recruitment = 1
    ),
    NA
  )
  expect_equal(
    plot_res_scaled$data$value[plot_res_scaled$data$stat_f == "Biomass"],
    plot_res$data$value[plot_res$data$stat_f == "Biomass"] * 1000
  )
  expect_equal(
    plot_res_scaled$data$value[plot_res_scaled$data$stat_f == "SSB"],
    plot_res$data$value[plot_res$data$stat_f == "SSB"] * 1000
  )
  expect_equal(
    plot_res_scaled$data$value[plot_res_scaled$data$stat_f == "Recruitment"],
    plot_res$data$value[plot_res$data$stat_f == "Recruitment"] * 1000
  )

  out_file <- tempfile("sam_ex_out_")
  expect_error(out_sam(res_test, filename = out_file), NA)
  expect_true(file.exists(paste0(out_file, ".csv")))

  par_file <- tempfile("sam_ex_par_")
  expect_error(par_table <- out_par(res_test, filename = par_file), NA)
  expect_true(file.exists(paste0(par_file, ".csv")))
  expect_true(is.data.frame(par_table))
  expect_true(all(c("FE", "MLE", "SE", "Gradient") %in% names(par_table)))
  expect_true(any(c("Unlinked value", "Unlinked.value") %in% names(par_table)))

  expect_error(index_plot2_res <- index_plot2(res_test), NA)
  expect_true(is.list(index_plot2_res))
  expect_true(all(c("index", "resid", "abund") %in% names(index_plot2_res)))

  expect_error(caa_plot_res <- caa_plot(res_test), NA)
  expect_true(is.list(caa_plot_res))
  expect_true(all(c("caa", "resid") %in% names(caa_plot_res)))

  expect_error(assess_res <- make_assess_result(res_test), NA)
  expect_true(is.data.frame(assess_res))
  expect_true(all(c("stat0", "Value", "SD", "CV", "lower", "upper", "Year", "Age", "stat", "Model") %in% names(assess_res)))
  expect_true(all(c("SSB", "Biomass", "F", "Catch", "Recruitment") %in% unique(assess_res$stat)))
})
