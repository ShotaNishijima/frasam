
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
  obj = samres$obj
  obs = samres$data$obs #fitするデータ (catch at age + index)
  if(is.null(subset)) subset <- 1:nrow(obs)
  osa.simple <- TMB::oneStepPredict(obj,
                                    observation.name = "logobs",
                                    method=method, #
                                    data.term.indicator = "keep",
                                    subset=1:nrow(obs),
                                    trace=trace)

  osa_resid = obs %>% as.data.frame() %>%
    mutate(logobs = log(obs)) %>%
    bind_cols(as.data.frame(osa.simple))

  return(osa_resid)
}

