source("R/sam.R")
source("R/boo_sam.R")
source("R/popsim_vpasam.R")
test_that <- function(desc, code) {
  force(code)
  cat(desc, ": OK\n")
}
expect_identical <- function(x, y) stopifnot(identical(x, y))
expect_equal <- function(x, y) stopifnot(isTRUE(all.equal(x, y)))
expect_true <- function(x) stopifnot(isTRUE(x))
source("tests/testthat/test-boo-sam-plusgroup.R")
source("tests/testthat/test-popsim-plusgroup.R")
