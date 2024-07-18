#' Popsimと同じ手法でSAMとVPAの疑似データを生成する
#'
#' @export


popsim_vpasam = function(Res,n=5,seed=1){
  if(is.numeric(set.seed)) set.seed(set.seed)
  
  if (class(Res)=="sam"){
    Res$input$bias.correct.sd <- FALSE
    Res0 <- Res
    dat.list <- list()
    pred.index <- Res$pred.index
    pred.caa <- Res$caa
    sigma.index = Res$sigma
    if (Res$input$est.method == "ls" && Res$input$index.key==NULL) sigma.index <- rep(Res$sigma,nrow(pred.index))
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