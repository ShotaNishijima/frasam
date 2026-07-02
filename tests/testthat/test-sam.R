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
  use_sam_tmb(overwrite=FALSE)
  args_def = formals(sam)
  input$cpp.file.name <- args_def$cpp.file.name
  input$p0.list <- NULL
  input$no_est <- TRUE
  testres <- safe_do_call(sam,input)

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
  use_sam_tmb(overwrite=FALSE)
  args_def = formals(sam)
  input$cpp.file.name <- args_def$cpp.file.name
  testres = safe_do_call(sam,input)
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
