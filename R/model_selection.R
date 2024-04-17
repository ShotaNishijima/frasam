
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
    minAIC2 = min(sapply(stage_sigma,function(X) X$aic))

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
#'
select_sigma_grid = function(
    samres,
    grid=expand.grid(var=c("varC","varF","varN","index.key")[1:2],X=1:2),
    stopAIC=TRUE){
  grid = grid %>% mutate(id = 1:n())

  samres2 = samres
  bestres = samres
  grid2 = grid
  minAIC = samres$aic
  reslist = list()
  # browser()
  # i<-2
  tbl_sigma = tibble("stage" = 0, "Age" = NA,"AIC" = samres2$aic)
  for(i in 1:nrow(grid)) {
    stage_sigma = lapply(1:nrow(grid2), function(j) {
      divide_sigma(samres2,var=as.character(grid2$var[j]),as.numeric(grid2$X[j]))
    })
    reslist[[i]] <- stage_sigma
    age_stage = which.min(sapply(stage_sigma,function(X) X$aic))
    tbl_sigma = bind_rows(tbl_sigma,tibble("stage" = i, "var" = as.character(grid2$var),"which" = as.numeric(grid2$X),"AIC" = sapply(stage_sigma,function(X) X$aic)))
    minAIC2 = min(sapply(stage_sigma,function(X) X$aic))

    message(paste0("Stage ", i, ": The selected setting is 'var='",as.character(grid2$var[age_stage]), " and 'which'=",as.numeric(grid2$X[age_stage]), " AIC = ", round(samres2$aic,2)))
    grid2 <- grid2[-age_stage,]
    samres2 = stage_sigma[[age_stage]]
    if (minAIC2 < minAIC) {
      bestres <- samres2
    }
    if (isTRUE(stopAIC) && minAIC2 >= minAIC) {
      # message(paste0("In the best setting, AIC = ", round(bestres$aic,2)))
      # message(paste0("The best setting is (", paste0(bestres$input[[var]],collapse=", "), "). AIC = ", round(bestres$aic,2),"\n"))
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
  # browser()
  message(paste0("In the best setting, AIC = ",round(min(tbl_sigma$AIC),2)))

  return(list(tbl_sigma=tbl_sigma,bestres = bestres,reslist=reslist))
}

