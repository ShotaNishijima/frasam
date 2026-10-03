library(frasam)
# devtools::load_all()

context("retro check")
test_that("test retrot",{
  # samres = get(load(system.file("data","samres_example.rda",package="frasam")))
  data("samres_example",package="frasam")
  # retrores = get(load(system.file("tests/testthat/testdata","retrores_example.rda",package="frasam")))
  retrores <- get(load(testthat::test_path("testdata", "retrores_example.rda")))
  input = samres$input
  input$p0.list <- NULL
  ensure_sam_tmb_loaded()
  args_def = formals(sam)
  input$cpp.file.name <- args_def$cpp.file.name
  testres = safe_do_call(sam,input)
  testretro = retro_sam(testres,n=1)

  dat_overwrite <- list(testretro$Res[[1]]$input$dat)
  idx <- which(!is.na(dat_overwrite[[1]]$index), arr.ind = TRUE)[1, ]
  dat_overwrite[[1]]$index[idx[1], idx[2]] <- dat_overwrite[[1]]$index[idx[1], idx[2]] * 1.01
  testretro_overwrite <- retro_sam(testres, n = 1, dat_overwrite = dat_overwrite)
  expect_equal(testretro_overwrite$Res[[1]]$input$dat, dat_overwrite[[1]])

  bad_dat_overwrite <- dat_overwrite
  bad_dat_overwrite[[1]]$index <- bad_dat_overwrite[[1]]$index[, -1, drop = FALSE]
  expect_error(
    retro_sam(testres, n = 1, dat_overwrite = bad_dat_overwrite),
    "same dimensions as the peeled data",
    fixed = TRUE
  )

  testcontents <-c("n","b","s","r","f")
  for(i in 1:length(testcontents)){
    expect_equal(eval(parse(text=paste0("retrores$retro.",testcontents[i],"[1]"))),
                 eval(parse(text=paste0("testretro$retro.",testcontents[i]))),tolerance = 1e-3)
  }
  for(i in 1:length(testcontents)){ #forecast
    expect_equal(eval(parse(text=paste0("retrores$retro.",testcontents[i],"2[1]"))),
                 eval(parse(text=paste0("testretro$retro.",testcontents[i],"2"))),tolerance = 1e-3)
  }
})
