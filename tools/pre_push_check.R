# Run before committing or pushing.
#
# Default:
#   Rscript tools/pre_push_check.R
#
# Include slow integration tests:
#   Rscript tools/pre_push_check.R --slow
# or:
#   Sys.setenv(FRASAM_RUN_SLOW_TESTS = "true")
#   source("tools/pre_push_check.R")

args <- commandArgs(trailingOnly = TRUE)
run_slow_tests <- "--slow" %in% args ||
  identical(Sys.getenv("FRASAM_RUN_SLOW_TESTS"), "true")

if (run_slow_tests) {
  Sys.setenv(FRASAM_RUN_SLOW_TESTS = "true")
  message("Running checks with slow integration tests enabled.")
} else {
  Sys.unsetenv("FRASAM_RUN_SLOW_TESTS")
  message("Running checks with slow integration tests skipped.")
  message("Use `Rscript tools/pre_push_check.R --slow` to include them.")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required. Install it with install.packages('devtools').")
}

devtools::load_all()
devtools::document()
devtools::test()
# devtools::test(filter = "sam")

devtools::check(vignettes = FALSE)
#vignetteの作成をスキップする場合

devtools::check()


## check building vignettes ----

devtools::load_all()

rmarkdown::render("vignettes/sam.Rmd", output_format = "rmarkdown::html_vignette")
rmarkdown::render("vignettes/FAQ.Rmd", output_format = "rmarkdown::html_vignette")

# pkgdown サイト全体を GitHub Actions に近い形で確認するなら:
pkgdown::build_site(new_process = FALSE, install = FALSE)

# パッケージの vignette として確認するなら:

devtools::build_vignettes()

