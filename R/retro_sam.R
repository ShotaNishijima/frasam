
#' レトロスペクティブ解析を実施する
#'
#' retrospective forecasting is also possible
#'
#' @param res SAM object
#' @param n the number of peels
#'
#' @export

retro_sam <- function(res, n=5, stat="mean", b.fix=TRUE,remove_short_index=-1, map_add = NULL, p0_retro_list = NULL){
  res.c <- res
  res.c$input$bias.correct.sd = FALSE
  Res <- list()
  obj.n <- obj.b <- obj.s <- obj.r <- obj.f <- NULL
  obj.n2 <- obj.b2 <- obj.s2 <- obj.r2 <- obj.f2 <- NULL
  max.a <- nrow(res$naa)

  if ("rec_logb" %in% names(map_add)) {
    res.c$input$b.init <- as.numeric(res$rec.par["b"])
  }

  if (isTRUE(b.fix)){
    # res.c$input$b.fix <- res$b
    res.c$input$b.fix <- exp(res$obj$env$parList()[["logB"]])
    # res.c$input$b.est <- FALSE
  }

  # i <- 1
  for (i in 1:n){
    nc <- ncol(res.c$input$dat$caa)

    res.c$input$dat$caa <- res.c$input$dat$caa[,-nc]
    for(j in 1:nrow(res.c$input$dat$index)){
      if (is.na(res.c$input$dat$index[j,ncol(res.c$input$dat$index)])) {
        res.c$input$dat$index[j,ncol(res.c$input$dat$index)-1] <- NA
      }
    }
    res.c$input$dat$index <- res.c$input$dat$index[,-nc,drop=FALSE]

    res.c$input$dat$catch.prop <- res.c$input$dat$catch.prop[,-nc]
    # res.c$input$catch_prop <- res.c$input$catch_prop[,-nc,]

    nc2 <- nc
    if (isTRUE(res$input$last.catch.zero)){
      res.c$input$dat$caa[,ncol(res.c$input$dat$caa)] <- 0
      nc2 <- nc2-1
    }

    res.c$input$p0.list <- res.c$par_list
    use.index = 1:nrow(res.c$input$dat$index)
    # res.c$input$retro.years <- i
    if (remove_short_index>0) {
      index_n = apply(res.c$input$dat$index,1,function(x) length(x)-sum(is.na(x)))
      # use.index = 1:nrow(res.c$input$dat$index)
      if (is.null(res.c$input$use.index)) {
        use.index = use.index[index_n > remove_short_index]
      } else {
        use.index = intersect(res.c$input$use.index,use.index[index_n > remove_short_index])
      }
      res.c$input$use.index <- use.index
      if (!is.null(res.c$input$index.key)) {
        res.c$input$index.key <- res.c$input$index.key[use.index]-min(res.c$input$index.key[use.index])+1
      }
      if (!is.null(res.c$input$index.b.key)) {
        res.c$input$index.b.key <- res.c$input$index.b.key[use.index]-min(res.c$input$index.b.key[use.index])
      }
      res.c$input$b.fix <- res.c$input$b.fix[unique(res.c$input$index.b.key)+1]
      if (!is.null(res.c$input$p0.list)) {
        res.c$input$p0.list$logB <- res.c$input$p0.list$logB[unique(res.c$input$index.b.key)+1]
      }
      res.c$input$p0.list$logSdLogObs <- c(res.c$input$p0.list$logSdLogObs[unique(res.c$input$varC)+1],
                                           res.c$input$p0.list$logSdLogObs[max(unique(res.c$input$varC))+1+unique(res.c$input$index.key)])

    }

    if (!is.null(map_add)) res.c$input$map.add <- map_add

    if (isTRUE(res.c$input$model_wm[1])) {
      if (i==1) res.c$input$weight_weight <- res.c$data$weight_weight
      res.c$input$weight_weight[,nc2] <- 0
    }
    if (isTRUE(res.c$input$model_wm[2])) {
      if (i==1) res.c$input$maturity_weight <- res.c$data$maturity_weight
      res.c$input$maturity_weight[,nc2] <- 0
    }

    if (is.null(p0_retro_list)) {
      res1 <- try(do.call(sam,res.c$input),silent=TRUE)
    } else {
      res.c$input$p0.list <- p0_retro_list[[i]]
      res1 <- try(do.call(sam,res.c$input),silent=TRUE)
    }

    if (class(res1) == "try-error") {
      res.c$input$p0.list <- res.c$par_list
      res.c$input$p0.list[["logQ"]] <- res.c$par_list[["logQ"]][use.index]
      res.c$input$p0.list[["logB"]] <- res.c$par_list[["logB"]][use.index]
      # res.c$input$b.fix %>% log
      res.c$input$p0.list[["logSdLogObs"]] <- c(unique(log(res.c$sigma.logC)),unique(log(res.c$sigma)[use.index]))
      res1 <- try(do.call(sam,res.c$input),silent=TRUE)
    }
    if (class(res1) == "try-error") {
      if (i>1) {
        res.c$input$p0.list <- res2$par_list
        res.c$input$p0.list[["logQ"]] <- res2$par_list[["logQ"]][use.index]
        res.c$input$p0.list[["logB"]] <- res2$par_list[["logB"]][use.index]
        # res.c$input$b.fix %>% log
        res.c$input$p0.list[["logSdLogObs"]] <- c(unique(log(res2$sigma.logC)),unique(log(res2$sigma)[use.index]))
        res1 <- try(do.call(sam,res.c$input),silent=TRUE)
      }
    }
    if (class(res1) == "try-error") {
      res.c$input$p0.list <- NULL
      res1 <- try(do.call(sam,res.c$input))
      if (class(res1) == "try-error") {
        stop(paste0("Error at trial #",i))
      }
    }

    Res[[i]] <- res1
    res2 <- res1
    if (res$input$last.catch.zero) Y <- nc-2 else Y <- nc-1

    obj.n <- c(obj.n, (sum(res1$naa[,Y])-sum(res$naa[,Y]))/sum(res$naa[,Y]))
    obj.b <- c(obj.b, (sum(res1$baa[,Y])-sum(res$baa[,Y]))/sum(res$baa[,Y]))
    obj.s <- c(obj.s, (sum(res1$ssb[,Y])-sum(res$ssb[,Y]))/sum(res$ssb[,Y]))
    obj.r <- c(obj.r, (res1$naa[1,Y]-res$naa[1,Y])/res$naa[1,Y])
    obj.f <- c(obj.f, (sum(res1$faa[-max.a,Y])-sum(res$faa[-max.a,Y]))/sum(res$faa[-max.a,Y]))

    # retrospective forecasting
    obj.n2 <- c(obj.n2, (sum(res1$naa[,Y+1])-sum(res$naa[,Y+1]))/sum(res$naa[,Y+1]))
    obj.b2 <- c(obj.b2, (sum(res1$baa[,Y+1])-sum(res$baa[,Y+1]))/sum(res$baa[,Y+1]))
    obj.s2 <- c(obj.s2, (sum(res1$ssb[,Y+1])-sum(res$ssb[,Y+1]))/sum(res$ssb[,Y+1]))
    obj.r2 <- c(obj.r2, (res1$naa[1,Y+1]-res$naa[1,Y+1])/res$naa[1,Y+1])
    obj.f2 <- c(obj.f2, (sum(res1$faa[-max.a,Y+1])-sum(res$faa[-max.a,Y+1]))/sum(res$faa[-max.a,Y+1]))
  }

  mohn <- c(get(stat)(obj.n,na.rm=TRUE),get(stat)(obj.b,na.rm=TRUE),get(stat)(obj.s,na.rm=TRUE),get(stat)(obj.r,na.rm=TRUE),get(stat)(obj.f,na.rm=TRUE))
  mohn2 <- c(get(stat)(obj.n2,na.rm=TRUE),get(stat)(obj.b2,na.rm=TRUE),get(stat)(obj.s2,na.rm=TRUE),get(stat)(obj.r2,na.rm=TRUE),get(stat)(obj.f2,na.rm=TRUE))

  names(mohn) <- names(mohn2) <- c("N","B","SSB","R","F")

  return(list(Res=Res,retro.n=obj.n, retro.b=obj.b, retro.s=obj.s, retro.r=obj.r, retro.f=obj.f, mohn=mohn,
              retro.n2=obj.n2, retro.b2=obj.b2, retro.s2=obj.s2, retro.r2=obj.r2, retro.f2=obj.f2, mohn_forecast=mohn2))
}


#' レトロの結果から各Indexに対する予測値を抽出して、Mean Absolute Scaled Errorを計算する関数
#'
#'
#' @param samres SAM object
#' @param retrores \code{retro_sam(res,...)}で実地されたレトロ解析の結果オブジェクト
#' @param log MASEを計算するときに、Indexの観測値と予測値に対してlogを取るかどうか（デフォルトはFALSE）
#' @param index_name 各Indexの名前ベクトル、Indexの数だけ必要
#'
#' @export

calc_mase = function(samres,
                     retrores,
                     h = 1,
                     log = FALSE,
                     index_name = NULL) {

  if (!is.null(index_name)) {
    if(length(index_name) != nrow(samres$input$dat$index)) stop("'length(index_name)' doesn't match with the number of indices")
  } else {
    index_name = str_c("Index ",as.character(1:nrow(samres$input$dat$index)))
  }

  convert_idx2tbl = function(x,value_name) {
    x %>% rownames_to_column(var="idx") %>%
      pivot_longer(cols = -idx, names_to = "year", values_to = value_name) %>%
      mutate(idx = as.integer(idx), year = as.integer(year)) %>%
      mutate(index = index_name[idx]) %>%
      mutate(index = fct_inorder(index))
  }

  obsdat = convert_idx2tbl(samres$input$dat$index,value_name="obs") %>%
    na.omit()

  fullpred = convert_idx2tbl(samres$pred.index,value_name="pred_full")

  basedat = left_join(obsdat,fullpred) %>% suppressMessages()
  nretro = length(retrores$Res)

  # make basic data of T and h for each index
  info = obsdat %>% group_by(idx, index) %>%
    summarise(T_i = max(year),
              T0_i = min(year),
              h = h)
  #

  cv_grid = expand.grid(idx = info$idx, retro_id = 1:nretro) %>%
    left_join(info) %>%
    mutate(
      year_target = T_i - retro_id + h,
      year_cond = T_i - retro_id
    )

  cv_grid = left_join(
    cv_grid,
    basedat %>% rename(year_target = year)
  ) %>% left_join(
    obsdat %>% rename(year_cond = year, obs_cond = obs)
  )

  cv_res <- dat_removed <- data.frame()

  for(i in 1:nretro) {
    res2 = retrores$Res[[i]]
    conddat = convert_idx2tbl(res2$pred.index,value_name="pred_cond")

    cv_res = bind_rows(
      cv_res,
      cv_grid %>%
      filter(retro_id==i) %>%
      left_join(conddat  %>% rename(year_target = year))
    ) %>% suppressMessages()

    conddat2 = conddat %>%
      left_join( cv_grid %>% filter(retro_id==i) %>% select(idx,year_target)) %>%
      filter(year <= year_target) %>%
      mutate(retro_id = i) %>% select(retro_id,everything()) %>%
      suppressMessages()

    dat_removed = bind_rows(
      dat_removed,
      conddat2
    )
  }

  if(isTRUE(log)) {
    cv_res = cv_res %>%
      mutate(obs = log(obs), obs_cond = log(obs_cond),
             pred_full = log(pred_full), pred_cond = log(pred_cond))
  }
  cv_res = cv_res %>%
    mutate(error_denom = obs - obs_cond,
           error_numer = obs - pred_cond) %>%
    select(retro_id, everything())

  mase_res = cv_res %>% group_by(idx, index) %>%
    summarise(denominator = mean(abs(error_denom), na.rm = TRUE),
              numerator = mean(abs(error_numer), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(MASE = numerator / denominator)

  return(list(full = basedat,
              removed = dat_removed,
              cv = cv_res,
              mase = mase_res))
}
