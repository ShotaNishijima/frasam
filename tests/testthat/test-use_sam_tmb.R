library(frasam)

context("use_sam_tmb check")

test_that("oututput value check",{
  res <- use_sam_tmb()
  expect_equal(res,TRUE)
})
