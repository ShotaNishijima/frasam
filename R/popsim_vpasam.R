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
    catch.sim <- sam_catch_simulation_inputs(Res)
    pred.caa <- catch.sim$mean
    sigma.index = Res$sigma
    if(length(sigma.index)==1) {
      sigma.index <- rep(sigma.index, nrow(Res$input$dat$index))
    }
    sigma.caa <- catch.sim$sd
    for (j in 1:n) {
      sim.dat <- Res$input$dat
      for (i in 1:nrow(pred.index)) {
        sim.index <- exp(rnorm(ncol(pred.index),as.numeric(log(pred.index[i,])),sigma.index[i]))
        sim.dat$index[i,!is.na(sim.dat$index[i,])] <- sim.index[!is.na(sim.dat$index[i,])]
      }
      for (i in 1:nrow(pred.caa)) {
        present <- catch.sim$present[i, ]
        sim.dat$caa[i,present] <- exp(rnorm(sum(present),
          log(pred.caa[i,present]),sigma.caa[i,present]))
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
    if(length(sigma.index)==1) {
      sigma.index <- rep(sigma.index, nrow(Res$input$dat$index))
    }
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
#' @param reest_lambda ridgeVPAで毎回lambdaを再推定するかどうか．TRUEの場合はlambda探索を実行する
#' @inheritParams frasyr::autocalc_ridgevpa target_retro n_retro b_fix
#' @param ... frasyr::autocalc_ridgevpaで設定可能な引数
#'
#' @export
fit2PSdata = function(res,
                      PSdata,
                      reest_lambda = FALSE,
                      target_retro = NULL,
                      n_retro = 5,
                      b_fix = TRUE,
                      ...) {

  if (!isTRUE(reest_lambda) && ...length() > 0L) {
    warning("'...' is ignored when reest_lambda = FALSE")
  }

  if (isTRUE(reest_lambda)) {
    if (inherits(res, "vpa")) {
      missing_args <- c(
        if (is.null(target_retro)) "target_retro",
        if (is.null(n_retro)) "n_retro",
        if (is.null(b_fix)) "b_fix"
      )
      if (length(missing_args) > 0) {
        stop(
          "When 'reest_lambda = TRUE', the following arguments must be specified: ",
          paste(missing_args, collapse = ", ")
        )
      }
    }
    if (inherits(res, "sam")) {
        warning("The 'reest_lambda' option is only for vpa so will be ignored")
      }
  }

  if (isTRUE(reest_lambda) &&
      inherits(res, "vpa") &&
      !is.null(res$input$eta)) {
    stop(
      "fit2PSdata() does not support re-estimating lambda ",
      "when eta is specified.",
      call. = FALSE
    )
  }

  res_list <- lapply(1:length(PSdata), function(i){
    input <- res$input
    input$dat <- PSdata[[i]]
    if(class(res)[1]=="vpa") {
      if (reest_lambda) {
                input$lambda <- 0
                res1 <- try(
                    do.call(vpa, input)
                )
                if(class(res1)[1]=="try-error"){
                  res2 <- res1
                } else {
                  ridge_res <- try(
                    frasyr::autocalc_ridgevpa(
                      input = res1$input,
                      target_retro = target_retro,
                      n_retro = n_retro,
                      b_fix = b_fix,
                      ...
                    )
                  )
                  if (class(ridge_res)[1]=="try-error") {
                    res2 <- ridge_res
                  } else {
                    lambda1 <- ridge_res$min_penalty
                    input$lambda <- lambda1
                    res2 <- try(
                      do.call(vpa, input)
                    )
                  }
                }
      } else {
                res2 <- try(
                    do.call(vpa, input)
                )
            }
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



