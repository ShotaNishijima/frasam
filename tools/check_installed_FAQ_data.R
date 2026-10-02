# Verify package data from an installed copy in a fresh R process.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  lib <- tempfile("frasam-install-check-")
  dir.create(lib)
  r <- file.path(R.home("bin"), "R.exe")
  status <- system2(r, c("CMD", "INSTALL", "--no-docs", "--no-multiarch",
                         paste0("--library=", shQuote(lib)), "."))
  stopifnot(status == 0L)
  status <- system2(file.path(R.home("bin"), "Rscript.exe"),
                    c(shQuote("tools/check_installed_FAQ_data.R"), shQuote(lib)))
  stopifnot(status == 0L)
} else {
  lib <- args[1L]
  library(frasam, lib.loc = lib)
  stopifnot(normalizePath(find.package("frasam")) ==
              normalizePath(file.path(lib, "frasam")))
  # Prevent data() from finding the source-tree data directory.
  setwd(tempdir())
  loaded <- new.env()
  data("res_bh", package = "frasam", lib.loc = lib, envir = loaded)
  data("res_vpa", package = "frasam", lib.loc = lib, envir = loaded)
  stopifnot(inherits(loaded$res_bh, "sam"),
            inherits(loaded$res_vpa, "vpa"),
            identical(loaded$res_bh$input$SR, "BH"))
  for (name in c("res_bh", "res_vpa")) {
    object <- loaded[[name]]
    stopifnot(is.matrix(as.matrix(object$naa)),
              any(is.finite(as.matrix(object$naa))))
    cat(name, ": class=", class(object), ", naa=",
        paste(dim(object$naa), collapse = " x "), "\n", sep = "")
  }
  stopifnot(inherits(frasam::res_bh, "sam"),
            inherits(frasam::res_vpa, "vpa"))
  cat("Verified installed package: data() and frasam:: access\n")
  cat("Installed path: ", find.package("frasam"), "\n", sep = "")
}
