fit_pseudo_sam_example <- function() {
  data("dat_ex", package = "frasam")

  use_sam_tmb(TmbFile = "sam2", overwrite = FALSE)

  sam(
    dat_ex,
    last.catch.zero = TRUE,
    abund = c("N", "N", "N", "SSB", "B"),
    min.age = c(0, 0, 1, 0, 0),
    max.age = c(0, 0, 1, 6, 6),
    rec.age = 0,
    index.key = 0:4,
    b.est = FALSE,
    SR = "RW",
    varC = c(0, 0, 1, 1, 1, 2, 2),
    varF = c(0, 0, 1, 1, 1, 1, 1),
    varN = c(0, 1, 1, 1, 1, 1, 1),
    rho.mode = 3,
    bias.correct = FALSE,
    silent = TRUE
  )
}

skip_if_not_slow_tests <- function() {
  testthat::skip_if(
    Sys.getenv("FRASAM_RUN_SLOW_TESTS") != "true",
    "Set FRASAM_RUN_SLOW_TESTS=true to run slow integration tests."
  )
}
