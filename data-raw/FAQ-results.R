# Recreate the example assessment results using the code in vignettes/FAQ.Rmd.
# Run from the package root with Rscript data-raw/FAQ-results.R.
if ("--check" %in% commandArgs(trailingOnly = TRUE)) {
  loaded <- new.env()
  data(list = c("res_vpa", "res_bh"), envir = loaded)
  stopifnot(
    is.list(loaded$res_vpa), is.list(loaded$res_bh),
    identical(dim(loaded$res_vpa$naa), dim(loaded$res_bh$naa)),
    identical(loaded$res_bh$input$SR, "BH"),
    any(is.finite(as.matrix(loaded$res_vpa$naa))),
    !any(is.infinite(as.matrix(loaded$res_vpa$naa))),
    all(is.finite(as.matrix(loaded$res_bh$naa)))
  )
  message("Verified data(): res_vpa and res_bh")
  quit(status = 0)
}

pkgload::load_all(export_all = FALSE)

faq <- readLines("vignettes/FAQ.Rmd", encoding = "UTF-8")
chunk_code <- function(label) {
  start <- grep(paste0("^```\\{r ", label, "[,}]"), faq)
  stopifnot(length(start) == 1L)
  end <- start + which(faq[(start + 1L):length(faq)] == "```")[1L]
  parse(text = faq[(start + 1L):(end - 1L)])
}

# Compile in a temporary directory so generated TMB files are not packaged.
run_dir <- tempfile("faq-tmb-")
dir.create(run_dir)
use_sam_tmb(RunDir = run_dir)

# The first five expressions load dat_ex, fit RW, and refit with BH.
for (expr in chunk_code("use-do.call")[1:5]) eval(expr)
# The first three expressions prepare dat_ex2 and fit the tuned VPA.
for (expr in chunk_code("set-qinit")[1:3]) eval(expr)

save(res_vpa, file = "data/res_vpa.rda", compress = "xz", version = 2)
save(res_bh, file = "data/res_bh.rda", compress = "xz", version = 2)

# Verify the public data() interface and the serialized numerical results.
loaded <- new.env(parent = baseenv())
data(list = c("res_vpa", "res_bh"), envir = loaded)
stopifnot(
  identical(loaded$res_vpa$naa, res_vpa$naa),
  identical(loaded$res_bh$naa, res_bh$naa),
  identical(loaded$res_bh$input$SR, "BH"),
  any(is.finite(as.matrix(loaded$res_vpa$naa))),
  !any(is.infinite(as.matrix(loaded$res_vpa$naa))),
  all(is.finite(as.matrix(loaded$res_bh$naa)))
)
message("Verified data(): res_vpa and res_bh")
