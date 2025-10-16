library(frasam)
# devtools::load_all()

context("retro check")
test_that("test retrot",{
  # samres = get(load(system.file("data","samres_example.rda",package="frasam")))
  data("samres_example",package="frasam")
  retrores = get(load(system.file("tests/testthat/testdata","retrores_example.rda",package="frasam")))
  input = samres$input
  input$p0.list <- NULL
  use_sam_tmb(overwrite=FALSE)
  args_def = formals(sam)
  input$cpp.file.name <- args_def$cpp.file.name
  testres = safe_do_call(sam,input)
  testretro = retro_sam(testres,n=1)

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
