library(frasam)

context("slow pseudo SAM workflows")

test_that("slow pseudo SAM workflows run", {
  skip_if_not_slow_tests()

  expect_error(res_test <- fit_pseudo_sam_example(), NA)

  data("retro_ex", package = "frasam")

  expect_error(
    retro_test <- suppressWarnings(retro_sam(res_test, n = 3)),
    NA
  )
  expect_equal(length(retro_test$Res), length(retro_ex$Res))

  retrocontents <- c(
    "retro.n", "retro.b", "retro.s", "retro.r", "retro.f",
    "retro.n2", "retro.b2", "retro.s2", "retro.r2", "retro.f2",
    "mohn", "mohn_forecast"
  )

  for (i in seq_along(retrocontents)) {
    expect_equal(
      eval(parse(text = paste0("retro_ex$", retrocontents[i]))),
      eval(parse(text = paste0("retro_test$", retrocontents[i]))),
      tolerance = 1e-2
    )
  }

  expect_error(
    retro_plot_res <- retro_plot(
      res_test,
      retro_test,
      plot_mohn = TRUE
    ),
    NA
  )
  expect_false(is.null(retro_plot_res))

  expect_error(
    capture.output(
      jitter_res <- do_jitter(
        res_test,
        nsim = 3,
        seed = 1
      )
    ),
    NA
  )
  expect_true(is.list(jitter_res))
  expect_true(all(c("resdat", "reslist") %in% names(jitter_res)))
  expect_equal(nrow(jitter_res$resdat), 4)
  expect_equal(length(jitter_res$reslist), 3)

  expect_error(
    osa_res <- do_osa_resid(
      res_test,
      subset = 101:103,
      trace = 0
    ),
    NA
  )
  expect_true(is.data.frame(osa_res))
  expect_equal(nrow(osa_res), 3)

  expect_error(osa_plot_res <- plot_osa_resid(osa_res), NA)
  expect_true(is.list(osa_plot_res))
  expect_true(all(c("caa", "index", "qq") %in% names(osa_plot_res)))

  expect_error(
    hindcast_plot_res <- plot_hindcastCV(
      res_test,
      retro_test,
      show_mase = TRUE
    ),
    NA
  )
  expect_false(is.null(hindcast_plot_res))

  expect_error(
    loo_res <- suppressWarnings(do_loo_index(res_test)),
    NA
  )
  expect_true(is.list(loo_res))
  expect_equal(length(loo_res), nrow(res_test$input$dat$index))
  loo_ok <- vapply(loo_res, function(x) inherits(x, "sam"), logical(1))
  loo_failed <- vapply(loo_res, function(x) inherits(x, "try-error"), logical(1))
  expect_true(all(loo_ok | loo_failed))
  expect_true(any(loo_ok))

  expect_error(
    capture.output(
      boot_res <- boo_sam(
        res_test,
        n = 10,
        seed = 1,
        method = "p"
      )
    ),
    NA
  )
  expect_true(is.list(boot_res))
  expect_equal(length(boot_res), 10)
  expect_true(all(vapply(boot_res, function(x) inherits(x, "sam"), logical(1))))

  expect_error(
    boosam_plot_res <- plot_boosam(
      res_test,
      boot_res,
      CI = 0.95
    ),
    NA
  )
  expect_false(is.null(boosam_plot_res))

  expect_error(
    capture.output(
      profile_res <- samprofile(
        res_test,
        param_name = "logQ",
        which_param = 1,
        length = 3
      )
    ),
    NA
  )
  expect_true(is.list(profile_res))
  expect_true(all(c("obj_tbl", "par_tbl", "obj_list", "opt_list") %in% names(profile_res)))
  expect_equal(length(profile_res$obj_list), 3)
  expect_equal(length(profile_res$opt_list), 3)

  expect_error(
    mixed_res <- est_mixed(
      tmbdata = res_test$data,
      par_init = res_test$par_list,
      map = res_test$map,
      random = unique(names(res_test$rep$par.random)),
      cpp.file.name = res_test$input$cpp.file.name,
      silent = TRUE,
      bias_correct = FALSE
    ),
    NA
  )
  expect_true(is.list(mixed_res))
  expect_true(all(c("obj", "opt", "rep") %in% names(mixed_res)))
})
