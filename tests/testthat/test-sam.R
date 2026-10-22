library(frasam)
# load_all()
context("input check")
test_that("test input",{
  # dat = get(load(system.file("data","dat_example.rda",package="frasam")))
  data("dat_example",package="frasam")
  # samres = get(load(system.file("data","samres_example.rda",package="frasam")))
  data("samres_example",package="frasam")
  input = samres$input
  tmbdata = samres$data
  ensure_sam_tmb_loaded()
  args_def = formals(sam)
  input$cpp.file.name <- args_def$cpp.file.name
  input$p0.list <- NULL
  input$no_est <- TRUE
  testres <- safe_do_call(sam,input)

  # Catch ages are internal row positions, independent of their labels.
  input_age1 <- input
  rownames(input_age1$dat$caa) <- as.character(seq_len(nrow(input_age1$dat$caa)))
  testres_age1 <- safe_do_call(sam, input_age1)
  expect_equal(testres_age1$data$obs, testres$data$obs)
  catch_obs <- testres_age1$data$obs[testres_age1$data$obs[, "fleet"] == 1, ]
  expect_equal(sort(unique(catch_obs[, "age"])),
               seq_len(nrow(input_age1$dat$caa)) - 1L)
  expect_equal(catch_obs[, "maxage"], catch_obs[, "age"])
  expect_equal(rownames(testres_age1$input$dat$caa), rownames(input_age1$dat$caa))

  input_remove <- input
  input_remove$remove.abund <- matrix(
    0,
    nrow = length(input_remove$abund),
    ncol = ncol(input_remove$dat$caa),
    dimnames = list(NULL, colnames(input_remove$dat$caa))
  )
  input_remove$remove.abund[2, ] <- 1
  testres_remove <- safe_do_call(sam, input_remove)
  expect_identical(testres_remove$data$use_remove_abund, 1L)
  expect_equal(ncol(testres_remove$data$obs), ncol(testres$data$obs) + 1)
  expect_true("remove" %in% colnames(testres_remove$data$obs))
  expect_equal(
    testres_remove$data$obs[testres_remove$data$obs[, "fleet"] == 1, "remove"],
    rep(0, sum(testres_remove$data$obs[, "fleet"] == 1))
  )
  expect_equal(
    testres_remove$data$obs[testres_remove$data$obs[, "fleet"] == 3, "remove"],
    rep(1, sum(testres_remove$data$obs[, "fleet"] == 3))
  )

  bad_input <- input
  bad_input$remove.abund <- matrix(0, nrow = length(bad_input$abund), ncol = ncol(bad_input$dat$caa) - 1)
  expect_error(
    safe_do_call(sam, bad_input),
    "The dimension of 'remove.abund'",
    fixed = TRUE
  )

  input_fix_zero_f <- input
  input_fix_zero_f$fix_zero_f <- matrix(
    0,
    nrow = nrow(input_fix_zero_f$dat$waa),
    ncol = ncol(input_fix_zero_f$dat$waa),
    dimnames = dimnames(input_fix_zero_f$dat$waa)
  )
  input_fix_zero_f$fix_zero_f[1, 1] <- 1
  zero_caa_year <- as.numeric(colnames(input_fix_zero_f$dat$caa)[1])
  input_fix_zero_f$dat$caa[1, 1] <- 0
  testres_fix_zero_f <- safe_do_call(sam, input_fix_zero_f)
  expect_equal(testres_fix_zero_f$data$fix_zero_f, input_fix_zero_f$fix_zero_f)
  zero_caa_obs <- testres_fix_zero_f$data$obs[
    testres_fix_zero_f$data$obs[, "fleet"] == 1 &
      testres_fix_zero_f$data$obs[, "age"] == 0 &
      testres_fix_zero_f$data$obs[, "year"] == zero_caa_year,
    ,
    drop = FALSE
  ]
  expect_equal(nrow(zero_caa_obs), 0)

  bad_input <- input
  bad_input$fix_zero_f <- matrix(
    0,
    nrow = nrow(bad_input$dat$waa),
    ncol = ncol(bad_input$dat$waa) - 1
  )
  expect_error(
    safe_do_call(sam, bad_input),
    "The dimension of 'fix_zero_f'",
    fixed = TRUE
  )

  bad_input <- input
  bad_input$q.init <- rep(1, length(bad_input$abund) + 1)
  expect_error(safe_do_call(sam,bad_input), "'q.init' must have length")

  bad_input <- input
  bad_input$p0.list <- testres$init
  bad_input$p0.list$logQ <- c(bad_input$p0.list$logQ, 0)
  expect_error(safe_do_call(sam,bad_input), "'p0.list' does not match the current model parameter structure")

  testcontents <-c("SR","b.est","b.fix","varC","varF","varN")
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste0("input$",testcontents[i]))),eval(parse(text=paste("testres$input$",testcontents[i]))))
  }

  # testres$data %>% names
  testcontents <-c("nobs","noYears","iy","minAge","maxAge","nlogF","nlogF")
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste("tmbdata$",testcontents[i]))),eval(parse(text=paste("testres$data$",testcontents[i]))))
  }
})


context("output check")
test_that("test output",{
  # dat = get(load(system.file("data","dat_example.rda",package="frasam")))
  data("dat_example")
  # samres = get(load(system.file("data","samres_example.rda",package="frasam")))
  data("samres_example")
  input = samres$input
  input$p0.list <- NULL
  # tmbdata = samres$data
  ensure_sam_tmb_loaded()
  args_def = formals(sam)
  input$cpp.file.name <- args_def$cpp.file.name
  testres = safe_do_call(sam,input)

  input_age1 <- input
  rownames(input_age1$dat$caa) <- as.character(seq_len(nrow(input_age1$dat$caa)))
  testres_age1 <- safe_do_call(sam, input_age1)
  expect_equal(testres_age1$opt$convergence, testres$opt$convergence)
  expect_equal(testres_age1$opt$objective, testres$opt$objective, tolerance = 1e-8)
  expect_equal(testres_age1$opt$par, testres$opt$par, tolerance = 1e-8)
  expect_equal(testres_age1$naa, testres$naa, tolerance = 1e-8)
  expect_equal(testres_age1$faa, testres$faa, tolerance = 1e-8)

  expect_equal(rownames(testres$naa), as.character(testres$data$minAge:testres$data$maxAge))
  expect_equal(rownames(testres$faa), as.character(testres$data$minAge:testres$data$maxAge))
  expect_equal(rownames(testres$caa), as.character(testres$data$minAge:testres$data$maxAge))
  expect_equal(colnames(testres$naa), as.character(testres$data$years))
  expect_equal(colnames(testres$faa), as.character(testres$data$years))
  expect_equal(colnames(testres$caa), as.character(testres$data$years))

  testcontents <-c("loglik","aic","q","b","opt$par","sigma","sigma.logC","sigma.logFsta","rho","phi","F","faa","N","naa")
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste0("samres$",testcontents[i]))),eval(parse(text=paste0("testres$",testcontents[i]))),tolerance = 1e-3)
  }

  ## test the case with caa matrix
  # pull request #36 https://github.com/ShotaNishijima/frasam/pull/36
  # test-pseudo-samよりcaaがmatrixの場合でもうまく行くかのテストをこちらに移行
  expect_equal(testres$opt$convergence,0)
  testres$input$dat$caa <- as.matrix(testres$input$dat$caa)
  testres$input$p0.list <- testres$par_list
  expect_true(is.matrix(testres$input$dat$caa))
  testres2 <- do.call(sam, testres$input)
  expect_equal(testres2$opt$convergence, 0)
  expect_true(is.data.frame(testres2$input$dat$caa))

  ## b.fix check with p0.list
  input_bfix <- testres$input
  expect_false(is.null(input_bfix$p0.list))
  input_bfix$b.fix[1] <- 2
  testres_bfix <- safe_do_call(sam, input_bfix)
  expect_equal(exp(testres_bfix$par_list$logB)[1], 2, tolerance = 1.0e-3)
  expect_equal(testres_bfix$b[1], 2, tolerance = 1.0e-3)

  ## b.fix check without p0.list
  input_bfix2 <- testres$input
  input_bfix2$b.fix[1] <- 2
  input_bfix2$p0.list <- NULL
  testres_bfix2 <- safe_do_call(sam, input_bfix2)
  expect_equal(exp(testres_bfix2$par_list$logB)[1], 2, tolerance = 1.0e-3)
  expect_equal(testres_bfix2$b[1], 2, tolerance = 1.0e-3)

  ## varN.fix check with p0.list
  input_varNfix <- testres$input
  input_varNfix$p0.list <- testres$par_list
  input_varNfix$varN.fix <- rep(NA, length(testres$par_list$logSdLogN))
  input_varNfix$varN.fix[1] <- 0.2^2
  testres_varNfix <- safe_do_call(sam, input_varNfix)
  expect_equal(exp(testres_varNfix$par_list$logSdLogN)[1], 0.2, tolerance = 1.0e-3)
  expect_equal(testres_varNfix$sigma.logN[1], 0.2, tolerance = 1.0e-3)

  ## bias-corrected catch-at-age should be consistent with reported catch biomass
  input_bias_correct <- testres$input
  input_bias_correct$p0.list <- testres$par_list
  input_bias_correct$bias.correct <- TRUE
  testres_bias_correct <- safe_do_call(sam, input_bias_correct)
  catch_biomass_sdr <- summary(testres_bias_correct$rep)
  catch_biomass_sdr <- catch_biomass_sdr[
    rownames(catch_biomass_sdr) == "Catch_biomass",
    "Est. (bias.correct)"
  ]
  catch_biomass_caa <- colSums(testres_bias_correct$caa * testres_bias_correct$waa_est)
  expect_lt(
    max(abs(catch_biomass_caa - catch_biomass_sdr) / abs(catch_biomass_sdr)),
    1.0e-5
  )

  ## rec.age > 0 should run without errors
  input_recage <- testres$input
  input_recage$rec.age <- 1
  input_recage$p0.list <- testres$par_list
  expect_warning(
    expect_error(testres_recage <- safe_do_call(sam, input_recage), NA),
    "rec.age > 0"
  )
  expect_equal(testres_recage$data$recAge, 1)
  expect_equal(as.numeric(testres_recage$data$minAge), 0)

  ## abund check
  # SSBに比例しているか
  expect_equal(sd(log(testres$pred.index[3,]/colSums(testres$ssb))),0,tolerance = 1.0e-3)
  expect_equal(sd(log(testres$pred.index[4,]/colSums(testres$ssb))),0,tolerance = 1.0e-3)
  # N, Bに比例しているか
  input = testres$input
  input$abund[3] <- "N"
  input$abund[4] <- "B"
  input$p0.list <- testres$par_list
  testres2 = safe_do_call(sam,input)
  expect_equal(sd(log(testres2$pred.index[3,]/colSums(testres2$naa))),0,tolerance = 1.0e-3)
  expect_equal(sd(log(testres2$pred.index[4,]/colSums(testres2$baa))),0,tolerance = 1.0e-3)

  input_remove = input
  input_remove$p0.list <- testres2$par_list
  input_remove$remove.abund <- matrix(
    0,
    nrow = length(input_remove$abund),
    ncol = ncol(input_remove$dat$caa),
    dimnames = list(NULL, colnames(input_remove$dat$caa))
  )
  input_remove$remove.abund[3, ] <- 1
  input_remove$remove.abund[4, ] <- input_remove$scale
  testres_remove = safe_do_call(sam,input_remove)
  expected_index3 <- testres_remove$q[3] *
    pmax(colSums(testres_remove$naa) - input_remove$remove.abund[3, ], 1e-6)
  expected_index4 <- testres_remove$q[4] *
    pmax(colSums(testres_remove$baa) / input_remove$scale - input_remove$remove.abund[4, ] / input_remove$scale, 1e-6)
  expect_equal(unname(as.numeric(testres_remove$pred.index[3, ])), unname(expected_index3), tolerance = 1.0e-6)
  expect_equal(unname(as.numeric(testres_remove$pred.index[4, ])), unname(expected_index4), tolerance = 1.0e-6)
  # Bs, Bfに比例しているか
  input = testres$input
  input$abund[3] <- "Bs"
  input$abund[4] <- "Bf"
  input$p0.list <- testres$par_list
  input$catch_prop <- testres$data$catch_prop4index
  input$catch_prop[1:7,,5] <- (0:6)/10
  # input$catch_prop[,,5]
  testres3 = suppressWarnings(safe_do_call(sam,input))
  expect_equal(sd(log(testres3$pred.index[3,]/colSums(testres3$baa*testres3$saa))),0,tolerance = 1.0e-3)
  expect_equal(sd(log(testres3$pred.index[4,]/colSums(testres3$baa*testres3$obj$report()[["saa_f"]][,,5]))),0,tolerance = 1.0e-3)

  ## ProcError check
  input = testres$input
  input$p0.list <- testres$par_list
  input$varN.fix <- NULL
  expect_error(testres4 <- safe_do_call(sam,input),NA)
  expect_false(all(testres4$sigma.logN ==  testres$sigma.logN))
  expect_equal(sd(testres4$sigma.logN[-1]),0)

  ## index.key = rep(0, nindex) should report the common index sigma for all indices
  input = testres$input
  input$index.key <- rep(0, length(input$abund))
  input$p0.list <- NULL
  expect_error(testres_ls <- safe_do_call(sam,input),NA)
  expect_equal(length(testres_ls$sigma), length(input$abund))
  expect_equal(sd(testres_ls$sigma), 0, tolerance = 1.0e-8)
  expect_false(any(is.na(testres_ls$sigma)))
})
