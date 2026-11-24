
#' 観測誤差とプロセス誤差をどこかの年齢間で分けて推定し直す関数
#'
#' @param samres sam object
#' @param var 分ける誤差の指定. "varC"（年齢別漁獲尾数の観測誤差）,"varF"（Fのプロセス誤差）,"varN"（Nのプロセス誤差）,"index.key"（指標の観測誤差）のいずれかを選択．
#' @param which どこで区切りをいれるか. 1だと0歳と1歳以上で分け（加入が0歳の場合）、指標値の場合は1本目と2本目以降の観測誤差を分ける
#' @export
#'
divide_sigma = function(
    samres,
    var=c("varC","varF","varN","index.key")[1],
    which = (1:6)[1]) {
  input = samres$input
  if (is.null(input[[var]]) || length(input[[var]])==1) {
    stop(paste0("The length of ", var, " should be the number of age classes or indices (only for 'index key')"))
  }
  KEY = input[[var]]
  KEY[(which+1):length(KEY)] <- KEY[(which+1):length(KEY)] + 1
  input[[var]] <- KEY
  if (!is.null(input$p0.list)) {
    if (var =="varC") {
      input$p0.list$logSdLogObs <- c(input$p0.list$logSdLogObs[1],input$p0.list$logSdLogObs)
    } else {
      if (var=="index.key") {
        input$p0.list$logSdLogObs <- c(input$p0.list$logSdLogObs,rev(input$p0.list$logSdLogObs)[1])
      } else {
        if(var=="varF") {
          input$p0.list$logSdLogFsta <- c(input$p0.list$logSdLogFsta,mean(input$p0.list$logSdLogFsta))
        } else {
          if (var=="varN") {
            input$p0.list$logSdLogN <- c(input$p0.list$logSdLogN,mean(input$p0.list$logSdLogN))
            if (!is.null(input$varN.fix)) {
              if (length(input$varN.fix) != length(input$p0.list$logSdLogN)) {
                stop("cannot divide 'varN' when the 'varN.fix' option is used")
              }
            }
          } else {
            stop("'var' is not identified")
          }
        }
      }
    }
  }
  RES = do.call(sam,input)
  return(RES)
}

#' 観測誤差やプロセス誤差のステップ形式のモデル選択（一つの変数について）
#'
#' \code{devide_sigma}関数を順々に実行し、AIC規準で最適なモデルを探索
#'
#' @inheritParams divide_sigma
#' @param X 境目を入れる場所の候補
#' @param stopAIC AICが小さくならなかった時点で計算をやめるか（default: TRUE)
#'
#' @export
#'
select_sigma = function(
    samres,
    var=c("varC","varF","varN","index.key")[1],
    X=NULL,
    stopAIC=TRUE){
  if(!is.null(X)) {
    if (var=="varF") {
      X = 1:(length(samres$input[[var]])-2) #最高齢のFは全年と同じなのでA-1とAでは分けない
    } else {
      X = 1:(length(samres$input[[var]])-1)
    }
  }

  samres2 = samres
  bestres = samres
  X2 = X
  minAIC = samres$aic
  reslist = list()
  tbl_sigma = tibble("stage" = 0, "Age" = NA,"AIC" = samres2$aic)
  for(i in 1:length(X)) {
    # browser()
    stage_sigma = lapply(X2, function(x) {
      divide_sigma(samres2,var=var,x)
    })
    reslist[[i]] <- stage_sigma
    age_stage = which.min(sapply(stage_sigma,function(X) X$aic))
    tbl_sigma = bind_rows(tbl_sigma,tibble("stage" = i, "Age" = X2,"AIC" = sapply(stage_sigma,function(X) X$aic)))
    minAIC2 = min(sapply(stage_sigma,function(X) X$aic),na.rm=T)

    X2 <- X2[-age_stage]
    samres2 = stage_sigma[[age_stage]]

    message(paste0("Stage ", i, ": The selected setting is (",paste0(samres2$input[[var]],collapse=", "), "). AIC = ", round(samres2$aic,2),""))
    if (minAIC2 < minAIC) {
      bestres <- samres2
    }
    if (isTRUE(stopAIC) && minAIC2 >= minAIC) {
      message(paste0("The best setting is (", paste0(bestres$input[[var]],collapse=", "), "). AIC = ", round(bestres$aic,2),"\n"))
      break
    }
    if (minAIC2 < minAIC) {
      minAIC <- minAIC2
    }
  }

  tbl_sigma = tbl_sigma %>% group_by(stage) %>%
    mutate(model = ifelse(AIC == min(AIC),"selected",NA)) %>%
    ungroup() %>%
    mutate(model = ifelse(AIC == min(AIC),"best",model))

  return(list(tbl_sigma=tbl_sigma,bestres = bestres,reslist=reslist))
}

#' 観測誤差やプロセス誤差のステップ形式のモデル選択（複数の変数について）
#'
#' \code{devide_sigma}関数を順々に実行し、AIC規準で最適なモデルを探索
#'
#' @inheritParams divide_sigma
#' @inheritParams select_sigma
#' @param grid 'var'と'X'からなるdata.frame
#' @param check_converge TRUE（既定値）の場合、\code{check_fit_sam}を通ったモデルのみ選択する。
#' @param SEmax 固定効果の標準誤差の上限。\code{check_fit_sam}の\code{se_max}に渡す。
#' @details 各ステップで適合チェックを通り、AICが有限の候補から最小AICのモデルを選ぶ。
#'   候補がなければ探索を終了する。初期モデルも同じ基準で判定し、合格モデルが
#'   一つもない場合は\code{bestres}を\code{NULL}とする。最適化の収束コードは
#'   \code{check_fit_sam}の既定動作と同様、合否判定には使わない。
#' @return \code{tbl_sigma}（\code{fit_ok}に適合チェックの合否を記録）、
#'   \code{bestres}、各ステップの候補を保存する\code{reslist}からなるリスト。
#'
#' @export
#'
select_sigma_grid = function(
    samres,
    grid=expand.grid(var=c("varC","varF","varN","index.key")[1:2],X=1:2),
    stopAIC=TRUE,
    check_converge=TRUE,
    SEmax = 10){
  fit_summary <- function(fit) {
    fit_ok <- check_fit_sam(fit, se_max = SEmax, verbose = FALSE)$ok
    tibble(
      AIC = fit$aic,
      convergence = if (is.null(fit$opt$convergence)) NA_integer_ else fit$opt$convergence,
      pdHess = if (is.null(fit$rep$pdHess)) NA else fit$rep$pdHess,
      maxSE = suppressWarnings(max(sqrt(diag(fit$rep$cov.fixed)), na.rm = TRUE)),
      fit_ok = isTRUE(fit_ok)
    )
  }
  eligible <- function(tbl) {
    is.finite(tbl$AIC) & (!isTRUE(check_converge) | tbl$fit_ok)
  }

  samres2 <- samres
  grid2 <- grid
  reslist <- list()
  tbl_sigma <- bind_cols(tibble(stage = 0L, var = NA_character_, which = NA_real_),
                         fit_summary(samres))
  tbl_sigma$model <- NA_character_
  best_row <- if (eligible(tbl_sigma)) 1L else NA_integer_
  bestres <- if (!is.na(best_row)) samres else NULL
  minAIC <- if (!is.na(best_row)) samres$aic else Inf

  for (i in seq_len(nrow(grid))) {
    stage_sigma <- lapply(seq_len(nrow(grid2)), function(j) {
      divide_sigma(samres2, var = as.character(grid2$var[j]), as.numeric(grid2$X[j]))
    })
    reslist[[i]] <- stage_sigma
    stage_tbl <- bind_cols(
      tibble(stage = i, var = as.character(grid2$var), which = as.numeric(grid2$X)),
      bind_rows(lapply(stage_sigma, fit_summary))
    )
    stage_tbl$model <- NA_character_
    candidates <- which(eligible(stage_tbl))
    if (!length(candidates)) {
      tbl_sigma <- bind_rows(tbl_sigma, stage_tbl)
      warning("No eligible models passed the fit checks with a finite AIC; stopping model selection.",
              call. = FALSE)
      break
    }

    age_stage <- candidates[which.min(stage_tbl$AIC[candidates])]
    minAIC2 <- stage_tbl$AIC[age_stage]
    stage_tbl$model[age_stage] <- "selected"
    selected_row <- nrow(tbl_sigma) + age_stage
    tbl_sigma <- bind_rows(tbl_sigma, stage_tbl)
    samres2 <- stage_sigma[[age_stage]]
    message(paste0("Stage ", i, ": The selected setting is 'var='",
                   as.character(grid2$var[age_stage]), " and 'which'=",
                   as.numeric(grid2$X[age_stage]), " AIC = ", round(minAIC2, 2)))
    grid2 <- grid2[-age_stage, , drop = FALSE]

    if (minAIC2 < minAIC) {
      bestres <- samres2
      best_row <- selected_row
    }
    if (isTRUE(stopAIC) && minAIC2 >= minAIC) break
    minAIC <- min(minAIC, minAIC2)
  }

  if (!is.na(best_row)) {
    tbl_sigma$model[best_row] <- "best"
    message(paste0("In the best setting, AIC = ", round(bestres$aic, 2)))
  }
  tbl_sigma <- tbl_sigma %>%
    dplyr::select(stage, var, which, convergence, pdHess, AIC, maxSE, fit_ok, model)
  return(list(tbl_sigma = tbl_sigma, bestres = bestres, reslist = reslist))
}
