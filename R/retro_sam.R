
####################################################
### Retrospective Analysis including Forecasting ###
####################################################

#'レトロスペクティブ解析を実施する
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

    nc2 <- nc
    if (res$input$last.catch.zero){
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

    # res1 <- do.call(sam,res.c$input)
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

