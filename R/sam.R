#########################################################
##### Core function of State-space Assessment model #####
#########################################################

sam <- function(dat,
                last.catch.zero = FALSE,
                cpp.file.name = "jsam9",
                tmb.run = FALSE,
                abund = c("N","N","SSB","SSB"),
                min.age = 0,
                index.age = rep(0,length(abund)),
                index.max.age = index.age,
                plus.group = TRUE,
                alpha = 1.0,
                b.est = FALSE,
                b.fix = rep(NA,length(abund)),
                varC = 0,
                varN = 0,
                varF = 0,
                varN.fix = NULL,
                est.method = "ml",
                SR = "BH", #"RW": random walk, "RI": Ricker
                AR = 0,
                rho.mode = 2,
                # MA = 0,
                q.init = NULL,
                sdFsta.init = NULL,
                sdLogN.init = NULL,
                sdLogObs.init = NULL,
                rho.init = NULL,
                a.init = NULL,
                b.init = NULL,
                ref.year = 1:5,
                bias.correct = TRUE,
                bias.correct.sd = FALSE,
                get.random.vcov = FALSE,
                silent = FALSE,
                remove.Fprocess.year = NULL,
                RW.Forder= 0,
                map.add = NULL,
                p0.list = NULL
                # retro.years = 0,
                # gamma = 10
){

  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  caa <- dat$caa
  if (last.catch.zero) caa[ncol(caa)] <- NULL

  maa <- dat$maa
  waa <- dat$waa
  M <- dat$M
  index <- dat$index

  obs <- cbind(expand.grid(as.numeric(rownames(caa)),as.numeric(colnames(caa))),1,unlist(caa))[,c(2,3,1,4)]

  nindex <- nrow(index)
  for (i in 1:nindex){
    index1 <- index[i,]
    index2 <- cbind(i+1,as.numeric(names(index1)[!is.na(index1)]),index.age[i],c(index1[!is.na(index1)]))[,c(2,1,3,4)]
    obs <- rbind(as.matrix(obs), index2)
  }

  colnames(obs) <- c("year","fleet","age","obs")
  rownames(obs) <- NULL

  obs <- as.matrix(obs)

  if (tmb.run) {
    library(TMB)
    compile(paste(cpp.file.name, ".cpp", sep = ""))
    dyn.load(dynlib(cpp.file.name))
  }

  data <- list()

  data$obs <- obs
  data$noFleets <- max(obs[,2])
  data$fleetTypes <- data$sampleTimes <- numeric(nindex+1)
  for (i in 1:nindex) {
    if (is.null(abund[i])) data$fleetTypes[i+1] <- data$freetTypes[i]
    if (abund[i] == "B") data$fleetTypes[i+1] <- 2
    if (abund[i] == "SSB") data$fleetTypes[i+1] <- 3
    if (abund[i] == "N") data$fleetTypes[i+1] <- 4
  }

  data$noYears <- ncol(waa)
  data$years <- as.numeric(colnames(waa))
  # if (retro.years>0) {
  #   for(i in 1:retro.years) data$years <- c(data$years,max(data$years)+1)
  # }
  data$iy <- data$obs[,1]-min(data$years) #year since the first year
  data$nobs <- nrow(data$obs)
  # colnames(maa) <- colnames(waa) <- colnames(M) <- range(data$obs[,1])[1]:((range(data$obs[,1])[2])+retro.years)
  data$propMat <- data$propMat2 <- as.matrix(t(maa))
  data$stockMeanWeight <- data$catchMeanWeight <- as.matrix(t(waa))
  data$natMor <- as.matrix(t(M))
  data$minAge <- as.matrix(min.age)
  data$maxAge <- as.matrix(min.age+nrow(caa)-1)
  data$maxAgePlusGroup <- ifelse(isTRUE(plus.group), 1, 0)
  data$rhoMode <- rho.mode
  ncol1 <- as.numeric(data$maxAge-data$minAge+1)
  data$landFrac <- data$disMeanWeight <- data$landMeanWeight <- data$propF <- data$propM <- matrix(0, nrow=data$noYears, ncol=ncol1)

  basemat <- matrix(-1, ncol = ncol1, nrow = max(data$obs[,2]))
  data$keyLogFsta <- data$keyLogQ <- data$keyLogB <- data$keyVarObs <- data$keyVarF <- data$keyVarLogN <- basemat

  data$keyLogFsta[1,] <- c(0:(data$maxAge-1),(data$maxAge-1))
  for (i in 1:nindex) data$keyLogQ[i+1,index.age[i]+1] <- i-1
  if (b.est) for (i in 1:nindex) data$keyLogB[i+1,index.age[i]+1] <- i-1

  data$keyVarObs[1,] <- varC
  if (est.method == "ls") {
    for (i in 1:nindex) data$keyVarObs[i+1,index.age[i]+1] <- max(varC)+1
  }
  if (est.method == "ml") {
    for (i in 1:nindex) data$keyVarObs[i+1,index.age[i]+1] <- max(varC)+i
  }

  data$keyVarF[1,] <- varF
  data$keyVarLogN[1,] <- varN

  if (SR == "BH") SR.mode <- 2 #Beverton-Holt
  if (SR == "RI") SR.mode <- 1 #Ricker
  if (SR == "RW") SR.mode <- 0 #Random walk
  if (SR == "HS") SR.mode <- 3 #Hockey-stick

  data$stockRecruitmentModelCode <- matrix(SR.mode)
  data$nlogF = max(data$keyLogFsta)+1
  data$nlogN = as.numeric(data$maxAge-data$minAge)+1
  # data$AR <- AR
  data$alpha <- alpha
  # data$MA <- MA
  # data$gamma <- gamma

  data$Fprocess_weight = rep(1,ncol(waa))
  if (!is.null(remove.Fprocess.year)) {
    data$Fprocess_weight[which(as.numeric(colnames(dat$waa)) %in% remove.Fprocess.year)] <- 0
  }
  data$F_RW_order <- RW.Forder


  U_init = rbind(
    matrix(5,nrow=ncol1,ncol=data$noYears),
    matrix(-1,nrow=max(data$keyLogFsta)+1,ncol=data$noYears)
    )

  logSdLogN_init = if (is.null(sdLogN.init)) rep(0.356675,max(data$keyVarLogN)+1) else log(sdLogN.init)

  if (!is.null(varN.fix)) {
    map_logSdLogN = 0:max(data$keyVarLogN)
    for(i in 1:(max(data$keyVarLogN)+1)) {
      if(!is.na(varN.fix[i])){
        logSdLogN_init[i] <- 0.5*log(varN.fix[i])
        map_logSdLogN[i] <- NA
      }
    }
  }

  if (is.null(p0.list)){
    params <- list(
      logQ      = if (is.null(q.init)) rep(-1,nindex) else log(q.init),
      logB      = sapply(1:nindex, function(i) ifelse(is.na(b.fix[i]), 0, log(b.fix[i]))),
      logSdLogFsta = if (is.null(sdFsta.init)) rep(-0.693147,max(data$keyVarF)+1) else log(sdFsta.init),
      logSdLogN    = logSdLogN_init,
      logSdLogObs  = if (is.null(sdLogObs.init)) rep(-0.356675,max(data$keyVarObs)+1) else log(sdLogObs.init),
      rec_loga     = if (is.null(a.init)) -4 else log(a.init),
      rec_logb     = if (is.null(b.init)) {if (SR=="HS") 13.5 else -14} else {log(b.init)},
      logit_rho  = if (is.null(rho.init)) 0 else log(rho.init/(1-rho.init)),
      # logScale     = numeric(data$noScaledYears),
      # logScaleSSB  = if(any(data$fleetTypes %in% c(3,4))) {numeric(0)} else {numeric(0)},
      # logPowSSB    = if(any(data$fleetTypes == 4))        {numeric(0)} else {numeric(0)},
      # logSdSSB     = if(any(data$fleetTypes %in% c(3,4))) {numeric(0)} else {numeric(0)},
      U = U_init,
      trans_phi1 = 0
    )
  } else {
    params <- p0.list
  }

  map <- list()
  if (!isTRUE(b.est)) map$logB <- factor(rep(NA,nindex)) else {
    if (!is.null(b.fix)) {
      map$logB <- 0:max(data$keyLogB)
      for (i in 1:nindex) if (!is.na(b.fix[i])) map$logB[i] <- NA
      map$logB <- factor(map$logB)
    }
  }

  if (SR=="RW") {
    map$rec_loga <- map$rec_logb <- factor(NA)
    map$trans_phi1 <- factor(NA)
  } else {
    if (AR==0) map$trans_phi1 <- factor(NA)
    # if (AR==1) map$phi2 <- factor(NA)
  }

  if (rho.mode==0) map$logit_rho <- factor(NA)
  if (rho.mode==1) map$logit_rho <- factor(NA)

  if(!is.null(varN.fix)) map$logSdLogN <- factor(map_logSdLogN)

  if (!is.null(map.add)) {
    tmp = messir::make_named_list(map.add)
    for(i in 1:length(tmp$map.add)) {
      map[[names(tmp$map.add)[i]]] <- map.add[[i]]
    }
  }

  obj <- MakeADFun(data, params, map = map, random=c("U"), DLL=cpp.file.name,silent=silent)

  lower <- obj$par*0-Inf
  upper <- obj$par*0+Inf

  opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper)
  if (opt$convergence!=0) warning("May not converge")
  rep <- sdreport(obj,bias.correct = bias.correct,bias.correct.control = list(sd=bias.correct.sd), getReportCovariance=get.random.vcov)
  if (max(rep$gradient.fixed)>1e-2) warning("Large maximum gradient component")

  SR.name <- "RW"
  is.SR <- ifelse(SR.mode==0, FALSE, TRUE)
  par_list = obj$env$parList(opt$par)

  BioRefPt <- function(data, rep, opt, ref.year, is.SR){
    if (is.null(rep$unbiased)) {
      NF <- matrix(rep$par.random,nrow=max(data$keyLogFsta)+1 + data$maxAge-data$minAge+1,ncol=data$noYears)
      logN <- NF[1:data$nlogN,]
      logF <- NF[(data$nlogN+1):(data$nlogN+data$nlogF),]
      logF <- logF[data$keyLogFsta[1,]+1,]
      logF[nrow(logF),] <- log(data$alpha)+logF[nrow(logF),]

      colnames(logN) <- colnames(logF) <- data$years
      rownames(logN) <- rownames(logF) <- data$minAge:data$maxAge
      naa <- exp(logN)
      faa <- exp(logF)
    } else {
      # logN <- matrix(rep$unbiased$value[names(rep$unbiased$value)=="logN"],ncol=data$noYears)
      # logF <- matrix(rep$unbiased$value[names(rep$unbiased$value)=="logF"],ncol=data$noYears)
      naa <- matrix(rep$unbiased$value[names(rep$unbiased$value)=="exp_logN"],ncol=data$noYears)
      faa <- matrix(rep$unbiased$value[names(rep$unbiased$value)=="exp_logF"],ncol=data$noYears)
      faa <- faa[data$keyLogFsta[1,]+1,]
      faa[nrow(faa),] <- data$alpha*faa[nrow(faa),]
      colnames(naa) <- colnames(faa) <- data$years
      rownames(naa) <- rownames(faa) <- data$minAge:data$maxAge
      logN <- log(naa)
      logF <- log(faa)
    }

    ref.year1 <- data$noYears-ref.year+1

    # naa <- exp(logN)
    # faa <- exp(logF)
    baa <- naa*t(data$stockMeanWeight)
    ssb <- baa*t(data$propMat)
    saa <- sweep(faa,2,apply(faa,2,max),FUN="/")
    zaa <- faa+t(data$natMor)
    caa <- faa/zaa*naa*(1-exp(-zaa))

    if (is.SR){
      # SR

      SR.name <- ifelse(as.numeric(data$stockRecruitmentModelCode)>1,
                        ifelse(as.numeric(data$stockRecruitmentModelCode)>2,"HS","BH"),
                        "RI")

      a <- exp(rep$par.fixed[names(rep$par.fixed)=="rec_loga"])
      b <- exp(rep$par.fixed[names(rep$par.fixed)=="rec_logb"])

      # B0

      rN0 <- rbind(1,exp(-apply(data$natM,1,cumsum)))
      rN0[data$nlogN,] <- rN0[data$nlogN,]/(1-rN0[data$nlogN+1,])
      rN0 <- rN0[1:data$nlogN,]

      ssb0 <- rN0*t(data$stockMeanWeight)*t(data$propMat)

      if (SR.name=="BH") R0 <- (a*colSums(ssb0)-1)/(b*colSums(ssb0))
      if (SR.name=="RI") R0 <- log(a*colSums(ssb0))/(b*colSums(ssb0))
      if (SR.name=="HS") R0 <- sapply(1:ncol(ssb0), function(i){
        uniroot(function(x){-2*x/a+x*colSums(ssb0)[i]+sqrt(1/b^2+data$gamma^2/4)-sqrt((x*colSums(ssb0)[i]-1/b)^2+data$gamma^2/4)}, c(0.001,max(naa[1,])*10))$root
        })

      N0 <- sweep(rN0,2,R0,FUN="*")

      B0 <- N0*t(data$stockMeanWeight)

      SSB0 <- sweep(ssb0,2,R0,FUN="*")

      mR0 <- mean(R0[ref.year1])
      mN0 <- sum(rowMeans(N0[,ref.year1]))
      mB0 <- sum(rowMeans(B0[,ref.year1]))
      mSSB0 <- sum(rowMeans(SSB0[,ref.year1]))

      BRP0 <- c(mR0,mN0,mB0,mSSB0)
      names(BRP0) <- c("R0","N0","B0","SSB0")

      # steepness

      if (SR.name == "BH") h <- 0.2*(1+b*mSSB0)/(1+0.2*b*mSSB0)
      if (SR.name == "RI") h <- 0.2*exp(0.8*b*mSSB0)
      if (SR.name == "HS") h <- 1-b/mSSB0

      rec.par <- c(a,b,h)
      names(rec.par) <- c("a","b","h")

      # MSY

      FAA <- exp(logF)
      SAA <- sweep(FAA, 2, colSums(FAA), FUN="/")
      mM <- rowMeans(t(data$natM)[,ref.year1])
      mSAA <- rowMeans(SAA[,ref.year1])
      mWAA <- rowMeans(t(data$stockMeanWeight)[,ref.year1])
      mMAA <- rowMeans(t(data$propMat)[,ref.year1])

      Fr <- exp(10)

      msy.func <- function(p){
        Fr <- exp(p)

        rN <- c(1,exp(-cumsum(mM+Fr*mSAA)))
        rN[data$nlogN] <- rN[data$nlogN]/(1-rN[data$nlogN+1])
        rN <- rN[1:data$nlogN]

        ssb1 <- rN*mWAA*mMAA

        if (SR.name=="BH") R1 <- (a*sum(ssb1, na.rm=TRUE)-1)/(b*sum(ssb1,na.rm=TRUE))
        if (SR.name=="RI") R1 <- log(a*sum(ssb1, na.rm=TRUE))/(b*sum(ssb1, na.rm=TRUE))
        if (SR.name=="HS") {
          R1 <- try(uniroot(function(x){-2*x/a+x*sum(ssb1, na.rm=TRUE)+sqrt(1/b^2+data$gamma^2/4)-sqrt((x*sum(ssb1, na.rm=TRUE)-1/b)^2+data$gamma^2/4)}, c(0.001,max(naa[1,])*10))$root,silent=TRUE)
          if (class(R1)=="try-error") R1 <- 0
        }

        zz <- mM+Fr*mSAA
        catch <- sum(Fr*mSAA/zz*R1*rN*mWAA*(1-exp(-zz)))
        return(-catch)
      }

      msy.res <- optimize(msy.func,c(-10,10))
      Fmsy <- exp(msy.res$minimum)
      MSY <- -msy.res$objective

      rN <- c(1,exp(-cumsum(mM+Fmsy*mSAA)))
      rN[data$nlogN] <- rN[data$nlogN]/(1-rN[data$nlogN+1])
      rN <- rN[1:data$nlogN]

      ssb.msy <- rN*mWAA*mMAA

      if (SR.name=="BH") Rmsy <- (a*sum(ssb.msy)-1)/(b*sum(ssb.msy))
      if (SR.name=="RI") Rmsy <- log(a*sum(ssb.msy))/(b*sum(ssb.msy))
      if (SR.name=="HS") Rmsy <- uniroot(function(x){-2*x/a+x*sum(ssb.msy)+sqrt(1/b^2+data$gamma^2/4)-sqrt((x*sum(ssb.msy)-1/b)^2+data$gamma^2/4)}, c(0.001,max(naa[1,])*10))$root

      Nmsy <- sum(Rmsy*rN)
      Bmsy <- sum(Rmsy*rN*mWAA)
      SSBmsy <- sum(Rmsy*rN*mWAA*mMAA)

      BRPmsy <- c(Fmsy, Nmsy, Bmsy, SSBmsy, Rmsy, MSY)
      names(BRPmsy) <- c("Fmsy","Nmsy", "Bmsy", "SSBmsy", "Rmsy", "MSY")
    } else {
      SR.rec <- rec.par <- BRP0 <- BRPmsy <- NULL
    }
    loglik <- -opt$objective
    aic <- 2*opt$objective+2*length(opt$par)

    q1 <- exp(as.numeric(rep$par.fixed[names(rep$par.fixed) == "logQ"]))
    if(!isTRUE(b.est)) b1 <- rep(1,nindex) else {
      if(is.null(b.fix)) b1 <- exp(as.numeric(rep$par.fixed[names(rep$par.fixed) == "logB"]))
      if(!is.null(b.fix)) {
        b1 <- b.fix
        b1[is.na(b1)] <- exp(as.numeric(rep$par.fixed[names(rep$par.fixed) == "logB"]))
      }
    }
    sigma0 <- exp(as.numeric(rep$par.fixed[names(rep$par.fixed) == "logSdLogObs"]))
    sigma1 <- sigma0[1:(max(data$keyVarObs[1,])+1)]
    sigma1 <- sapply(1:ncol1, function(i) sigma1[data$keyVarObs[1,i]+1])
    sigma2 <- sigma0[(max(data$keyVarObs[1,])+2):length(sigma0)]
    sigma3 <- exp(as.numeric(par_list[["logSdLogFsta"]]))
    sigma3 <- sapply(1:ncol1, function(i) sigma3[data$keyVarF[1,i]+1])
    sigma4 <- exp(obj$env$parList(opt$par)$logSdLogN)
    # sigma4 <- exp(as.numeric(rep$par.fixed[names(rep$par.fixed) == "logSdLogN"]))
    sigma4 <- sapply(1:ncol1, function(i) sigma4[data$keyVarLogN[1,i]+1])
    rho1 <- ifelse(as.numeric(data$rhoMode) > 1, 1/(1+exp(-as.numeric(rep$par.fixed[names(rep$par.fixed) == "logit_rho"]))),as.numeric(data$rhoMode))
    phi<-numeric(1)
    if (AR>0) phi[1] = rep$value[names(rep$value)=="phi1"][1]
    pred.index <- data.frame(matrix(0,nrow=nrow(dat$index),ncol=ncol(dat$waa)))
    colnames(pred.index) <- colnames(dat$waa)
    for (i in 1:nrow(pred.index)){
      if (abund[i]=="N") pred.index[i,] <- as.numeric(naa[index.age[i]+1,])
      if (abund[i]=="B") pred.index[i,] <- as.numeric(baa[index.age[i]+1,])
      if (abund[i]=="SSB") pred.index[i,] <- as.numeric(colSums(ssb))
      if (is.na(b1[i])) {
        pred.index[i,] <- q1[i]*pred.index[i,]
        } else {
          pred.index[i,] <- q1[i]*pred.index[i,]^b1[i]
          }
    }

    output <- list(data=data, obj=obj, opt=opt, rep=rep, q=q1, b=b1, sigma=sigma2, sigma.logC=sigma1, sigma.logFsta=sigma3, sigma.logN=sigma4, rho=rho1, phi=phi, loglik=loglik,aic=aic,N=exp(logN),F=exp(logF),ref.year=rev(ref.year1),SR=SR.name,rec.par=rec.par, BRP0=BRP0,BRPmsy=BRPmsy,naa=naa,faa=faa,baa=baa,ssb=ssb,saa=saa,zaa=zaa,caa=caa,pred.index=pred.index)
    return(output)
  }

  brp1 <- BioRefPt(data, rep, opt, ref.year, is.SR)
  brp1$input <- arglist
  brp1$par_list = obj$env$parList(opt$par)
  return(brp1)
}




