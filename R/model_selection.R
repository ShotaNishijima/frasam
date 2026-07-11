
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
#'
#' @export
#'
select_sigma_grid = function(
    samres,
    grid=expand.grid(var=c("varC","varF","varN","index.key")[1:2],X=1:2),
    stopAIC=TRUE,
    check_converge=FALSE,
    SEmax = 10){
  grid = grid %>% mutate(id = 1:n())

  samres2 = samres
  bestres = samres
  grid2 = grid
  minAIC = samres$aic
  reslist = list()
  # i<-2
  tbl_sigma = tibble("stage" = 0, "Age" = NA,"AIC" = samres2$aic,
                     "convergence" = samres2$opt$convergence,
                     "pdHess" = samres2$rep$pdHess,
                     "maxSE" = max(sqrt(diag(samres2$rep$cov.fixed)),na.rm=TRUE))
  for(i in 1:nrow(grid)) {
    stage_sigma = lapply(1:nrow(grid2), function(j) {
      divide_sigma(samres2,var=as.character(grid2$var[j]),as.numeric(grid2$X[j]))
    })
    reslist[[i]] <- stage_sigma

    age_stage = which.min(sapply(stage_sigma,function(X) X$aic))
    minAIC2 = min(sapply(stage_sigma,function(X) X$aic),na.rm=T)

    if( isTRUE(check_converge)) {
      # 収束していないもの、Hessianが求まっていないものを除くような仕様を追加(2025/05/14)
      convergence_aic = tibble(
        convergence = sapply(stage_sigma,function(X) X$opt$convergence),
        pdHess = sapply(stage_sigma,function(X) X$rep$pdHess),
        aic = sapply(stage_sigma,function(X) X$aic),
        maxSE = sapply(stage_sigma,function(X) max(sqrt(diag(X$rep$cov.fixed)),na.rm=TRUE))
      ) %>% mutate(ID = 1:n())

      tmp = convergence_aic %>%
        filter(convergence==0 & pdHess==TRUE & maxSE <= SEmax) %>%
        filter(aic == min(aic,na.rm=T))

      if (nrow(tmp)==0) {
        warning("Any models fail to converge")
        minAIC2 <- minAIC
        age_stage= tmp %>% pull(ID)
      } else {
        minAIC2 = tmp %>% pull(aic)
      }
        }

    tbl_sigma = bind_rows(tbl_sigma,
                          tibble("stage" = i, "var" = as.character(grid2$var),"which" = as.numeric(grid2$X),"AIC" = sapply(stage_sigma,function(X) X$aic),
                                 "convergence" = sapply(stage_sigma,function(X) X$opt$convergence),
                                 "pdHess" = sapply(stage_sigma,function(X) X$rep$pdHess),
                                 "maxSE" = sapply(stage_sigma,function(X) max(sqrt(diag(X$rep$cov.fixed)),na.rm=TRUE))))

    message(paste0("Stage ", i, ": The selected setting is 'var='",as.character(grid2$var[age_stage]), " and 'which'=",as.numeric(grid2$X[age_stage]), " AIC = ", round(minAIC2,2)))
    if (nrow(tmp)>0) {
      grid2 <- grid2[-age_stage,]
      samres2 = stage_sigma[[age_stage]]
    }
    # samres2$rep
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

  if(isTRUE(check_converge)) {
    tbl_sigma = tbl_sigma %>%
      mutate(SE_ok = ifelse(maxSE <= SEmax,TRUE,FALSE)) %>%
      group_by(stage,convergence,pdHess,SE_ok) %>%
      mutate(model = ifelse(AIC == min(AIC) & convergence==0 & pdHess == TRUE & SE_ok == TRUE,
                            "selected",NA)) %>%
      ungroup() %>%
      group_by(convergence,pdHess,SE_ok) %>%
      mutate(model = ifelse(AIC == min(AIC) & convergence==0 & pdHess == TRUE &  SE_ok == TRUE,"best",model)) %>%
      ungroup() %>% dplyr::select(-SE_ok)
  } else {
    tbl_sigma = tbl_sigma %>% group_by(stage) %>%
      mutate(model = ifelse(AIC == min(AIC),"selected",NA)) %>%
      ungroup() %>%
      mutate(model = ifelse(AIC == min(AIC),"best",model))
  }

  tbl_sigma = tbl_sigma  %>%
    dplyr::select(-Age,stage,var,which,convergence,pdHess,AIC,maxSE,model)
  # browser()
  message(paste0("In the best setting, AIC = ",round(filter(tbl_sigma,model=="best") %>% pull(AIC),2)))

  return(list(tbl_sigma=tbl_sigma,bestres = bestres,reslist=reslist))
}

