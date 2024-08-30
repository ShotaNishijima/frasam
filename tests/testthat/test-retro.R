library(frasam)
# devtools::load_all()

context("retro check")
test_that("test retrot",{
  samres = get(load(system.file("data","samres_example.rda",package="frasam")))
  retrores = get(load(system.file("tests/testthat/testdata","retrores_example.rda",package="frasam")))
  input = samres$input
  input$p0.list <- NULL
  use_sam_tmb(overwrite=FALSE)
  args_def = formals(sam)
  input$cpp.file.name <- args_def$cpp.file.name
  testres = safe_do_call(sam,input)
  # testres$input$p0.list <- testres$par_list
  testretro = retro_sam(testres,n=1)
  testretro = retro_sam(samres,n=1)

  testcontents <-c("loglik","aic","q","b","opt$par","sigma","sigma.logC","sigma.logFsta","rho","phi","F","faa","N","naa")
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste0("samres$",testcontents[i]))),eval(parse(text=paste0("testres$",testcontents[i]))),tolerance = 1e-3)
  }

  ## b.fix check
  input = testres$input
  input$b.fix[1] <- 2
  testres = do.call(sam,input)
  expect_equal(exp(testres$par_list$logB)[1],2,tolerance = 1.0e-3)
  expect_equal(testres$b[1],2,tolerance = 1.0e-3)

  ## abund check
  # SSBに比例しているか
  expect_equal(sd(log(testres$pred.index[3,]/colSums(testres$ssb))),0,tolerance = 1.0e-3)
  expect_equal(sd(log(testres$pred.index[4,]/colSums(testres$ssb))),0,tolerance = 1.0e-3)
  # N, Bに比例しているか
  input = testres$input
  input$abund[3] <- "N"
  input$abund[4] <- "B"
  input$p0.list <- testres$par_list
  testres2 = do.call(sam,input)
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
  testres3 = suppressWarnings(do.call(sam,input))
  expect_equal(sd(log(testres3$pred.index[3,]/colSums(testres3$baa*testres3$saa))),0,tolerance = 1.0e-3)
  expect_equal(sd(log(testres3$pred.index[4,]/colSums(testres3$baa*testres3$obj$report()[["saa_f"]][,,5]))),0,tolerance = 1.0e-3)

  ## ProcError check
  input = testres$input
  input$p0.list <- testres$par_list
  input$varN.fix <- NULL
  expect_error(testres4 <- do.call(sam,input),NA)
  expect_false(all(testres4$sigma.logN ==  testres$sigma.logN))
  expect_equal(sd(testres4$sigma.logN[-1]),0)
})
