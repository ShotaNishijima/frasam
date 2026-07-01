
#' SAMのOSA residualを計算する関数
#'
#' @param samres sam object
#'
#' @inheritParams TMB::oneStepPredict
#'
#' @encoding UTF-8
#'
#' @export

do_osa_resid <- function(samres,
                      method="oneStepGaussianOffMode", #
                      subset=NULL,
                      trace=2) {

  obj <- samres$obj
  obs_all <- samres$data$obs
  if (is.null(subset)) subset <- seq_len(nrow(obs_all))

  osa.simple <- TMB::oneStepPredict(
    obj,
    observation.name = "logobs",
    method = method,
    data.term.indicator = "keep",
    subset = subset,
    trace = trace
  )

  osa_resid <- obs_all[subset, , drop = FALSE] %>%
    as.data.frame() %>%
    mutate(logobs = log(obs)) %>%
    bind_cols(as.data.frame(osa.simple))

  return(osa_resid)
}

