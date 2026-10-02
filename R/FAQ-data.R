#' Example assessment results from the FAQ
#'
#' Assessment results generated using the examples in \code{vignettes/FAQ.Rmd}
#' and the packaged \code{dat_ex} data.
#'
#' @format Lists containing assessment inputs, estimated population numbers,
#'   fishing mortality, biomass, spawning biomass, and model diagnostics.
#' @details \code{res_vpa} is the tuned VPA result from the \code{set-qinit}
#'   example, with the final-year recruitment indices set to missing.
#'   \code{res_bh} is the SAM result with Beverton-Holt recruitment from the
#'   \code{use-do.call} example. Saved TMB external pointers cannot be used in
#'   a new R session; refit with \code{do.call(sam, res_bh$input)} when a live
#'   TMB objective is required, after running \code{use_sam_tmb()}.
#' @examples
#' data("res_vpa", package = "frasam")
#' data("res_bh", package = "frasam")
#' dim(res_vpa$naa)
#' dim(res_bh$naa)
#' @name FAQ-results
NULL

#' @rdname FAQ-results
"res_vpa"

#' @rdname FAQ-results
"res_bh"
