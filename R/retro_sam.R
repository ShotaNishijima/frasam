
####################################################
### Retrospective Analysis including Forecasting ###
####################################################

#'レトロスペクティブ解析を実施する
#'
#' @export

retro_sam <- function(res, n=5, stat="mean", b.fix=TRUE,remove_short_index=TRUE, map_add = NULL){
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
    res.c$input$b.fix <- res$b
    # res.c$input$b.est <- FALSE
  }

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

    if (res$input$last.catch.zero){
      res.c$input$dat$caa[,ncol(res.c$input$dat$caa)] <- 0
    }
    
    res.c$input$p0.list <- res.c$par_list
    # res.c$input$retro.years <- i
    if (isTRUE(remove_short_index)) {
      index_n = apply(res.c$input$dat$index,1,function(x) length(x)-sum(is.na(x)))
      use.index = 1:nrow(res.c$input$dat$index)
      if (is.null(res.c$input$use.index)) {
        use.index = use.index[index_n > 2] 
      } else {
        use.index = intersect(res.c$input$use.index,use.index[index_n > 2])
      }
      res.c$input$use.index <- use.index
    }
    
    if (!is.null(map_add)) res.c$input$map.add <- map_add

    # res1 <- do.call(sam,res.c$input)
    res1 <- try(do.call(sam,res.c$input),silent=TRUE)
    if (class(res1) == "try-error") {
      if (i>1) {
        res.c$input$p0.list <- res2$par_list
        res1 <- try(do.call(sam,res.c$input),silent=TRUE)
      }
    }
    if (class(res1) == "try-error") {
      res.c$input$p0.list <- NULL
      res1 <- do.call(sam,res.c$input)
    }

    Res[[i]] <- res1
    res2 <- res1
    Y <- nc-1

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

