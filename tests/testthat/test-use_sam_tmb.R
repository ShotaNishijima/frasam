library(frasam)

context("use_sam_tmb check")
test_that("test use_sam_tmb",{
  res <- use_sam_tmb()
  expect_equal(res,TRUE)
})

test_that("use_sam_tmb accepts auto_update argument", {
  res <- use_sam_tmb(auto_update = FALSE)
  expect_true(res)
})
