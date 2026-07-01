#' Popsimと同じ手法でSAMとVPAの疑似データを生成する
#'
#' @export


popsim_vpasam = function(Res,n=5,seed=1){
  if (!is.null(seed)) set.seed(seed)

  if (class(Res)=="sam"){
    Res$input$bias.correct.sd <- FALSE
    Res0 <- Res
    dat.list <- list()
    pred.index <- Res$pred.index
    pred.caa <- Res$caa
    sigma.index = Res$sigma
    if (Res$input$est.method == "ls" && is.null(Res$input$index.key)) sigma.index <- rep(Res$sigma,nrow(pred.index))
    sigma.caa <- Res$sigma.logC
    resid.index <- log(as.matrix(Res$input$dat$index))-log(as.matrix(Res$pred.index))
    resid.caa <- log(Res$input$dat$caa)-log(Res$caa)
    if (Res$input$last.catch.zero) resid.caa[,ncol(resid.caa)] <- NULL
    for (j in 1:n) {
      sim.dat <- Res$input$dat
      for (i in 1:nrow(pred.index)) {
        sim.index <- exp(rnorm(ncol(pred.index),as.numeric(log(pred.index[i,])),sigma.index[i]))
        sim.dat$index[i,!is.na(sim.dat$index[i,])] <- sim.index[!is.na(sim.dat$index[i,])]
      }
      for (i in 1:nrow(pred.caa)) {
        sim.dat$caa[i,] <- exp(rnorm(ncol(pred.caa),log(pred.caa[i,]),sigma.caa[i]))
      }
      if (Res$input$last.catch.zero) sim.dat$caa[,ncol(pred.caa)] <- 0
      dat.list[[j]] <- sim.dat
    }
  }
  if (class(Res)=="vpa"){
    Res0 <- Res
    dat.list <- list()
    pred.index <- Res$pred.index
    sigma.index = Res$sigma
    if (Res$input$est.method == "ls") sigma.index <- rep(Res$sigma,nrow(pred.index))
    resid.index <- log(as.matrix(Res$input$dat$index))-log(as.matrix(Res$pred.index))
    for (j in 1:n) {
      sim.dat <- Res$input$dat
      for (i in 1:nrow(pred.index)) {
        sim.index <- exp(rnorm(ncol(pred.index),as.numeric(log(pred.index[i,])),sigma.index[i]))
        sim.dat$index[i,!is.na(sim.dat$index[i,])] <- sim.index[!is.na(sim.dat$index[i,])]
      }
      dat.list[[j]] <- sim.dat
    }
  }

  return(dat.list)
}


#' 生成された疑似データにVPA/SAMを推定し、Self-test/Cross-testを実行
#'
#' @param res フィットさせるSAM/VPAのオブジェクト
#' @param PSdata \code{popsim_vpa()}で生成された疑似データ
#'
#' @export
fit2PSdata = function(res, PSdata) {
  res_list <- lapply(1:length(PSdata), function(i){
    input <- res$input
    input$dat <- PSdata[[i]]
    if(class(res)[1]=="vpa") {
      res2 <- try(do.call(vpa, input))
    } else {
      if(class(res)[1]=="sam") {
        res2 <- try(do.call(sam, input))
      } else {
        stop("'class(res)' should be either 'vpa' or 'sam'")
      }
    }
    return( res2 )
  })
  return( res_list )
}


#' Self-test, Cross-testの結果をまとめるための関数
#'
#' @param res_true 真の推定値をもつSAMまたはVPAのオブジェクト
#' @param fit2PS \code{fit2PSdata()} で得られる、疑似データにフィットさせた結果オブジェクト
#'
#' @inheritParams calc_metrics
#'
#' @return
#' A list with the following components:
#' \itemize{
#'   \item \code{summary}: A data.frame containing performance metrics.
#'     See [calc_metrics()] for details on the definitions of RMSE, MAE, R2_model, etc.
#'   \item \code{all}: すべてのRunについての結果のデータフレーム
#' }
#'
#' @encoding UTF-8
#'
#' @export
#
sumup_popsim <- function(
    res_true,
    fit2PS,
    percent = FALSE,
    CI = 0.95,
    ...) {
  tbl_true = convert_sam_tibble(res_true) %>%
    dplyr::select(-sim) %>%
    rename(value_true = value, type_true = type)
  tbl_ps = map_dfr(1:length(fit2PS), function(i) {
    convert_sam_tibble(fit2PS[[i]]) %>% mutate(simID = i)
  }) %>% dplyr::select(-sim)

  tbl_all = left_join(tbl_ps,tbl_true) %>% suppressMessages() %>%
    dplyr::select(simID,stat,year,age,value,value_true,type,type_true,everything())

  tbl_summary = tbl_all %>% group_by(stat,year,age) %>%
    summarise(metrics = list(calc_metrics(value_true, value, percent = percent, ...)),
              value_true = mean(value_true),
              lower = quantile(value,probs=0.5*(1-CI))[1],
              upper = quantile(value,probs=0.5*(1+CI))[1]
              ) %>%
    unnest_wider(metrics) %>%
    ungroup() %>% suppressMessages()

  return( list(all = tbl_all, summary = tbl_summary) )
}



