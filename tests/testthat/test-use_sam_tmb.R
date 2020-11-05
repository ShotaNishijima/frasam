library(frasam)

context("use_sam_tmb check")
test_that("test use_sam_tmb",{
  res <- use_sam_tmb()
  expect_equal(res,TRUE)
})
