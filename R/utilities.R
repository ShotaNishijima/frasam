#'
#' utilities for extracting statistics
#'
#'
#' @param res_sam result of sam function
#' @param res_future result of future projection
#' @param waa_catch weight at age for calculating total catch weight. This matrix should be given because SAM doesn't use this biological parameters
#' @param last_year last year of stock assessment
#' @param waa_biom weight at age for calculating biomass weight. If this is not given, sam_res$input$dat$waa is automatically used
#' @param year_biol year range to derive biological parameters. Biological parameters are averaged by age during the given period
#' @param year_Fcur year range to derive F at age considred as F current. Fs at age are averaged by age during the given period
#' @param perSPR percent SPR (%) for calculating F%SPR
#'
#' @examples
#' \donotrun{
#' res_pm_all <- purrr::map_dfr(basecase_list,
#'                   function(x) get_pm(x, NULL, waa_catch=x$input$dat$waa), .id="id") %>%
#'                   pivot_wider(values_from=value, names_from=stat)
#' }
#'
#' @export
#' 
#' @encoding UTF-8
#' 


# get performance measures
get_pm <- function(res_sam,
                   res_future=NULL,
                   waa_catch=res_sam$input$dat$waa,                     
                   last_year=2022,
                   waa_biom=res_sam$input$dat$waa,
                   year_biol=2020:2022,
                   year_Fcur=2020:2022,
                   perSPR=c(30,40,50, 60, 70)
                   ){

  # define year range
  last5year <- (last_year-(4:0)) %>% as.character()
  last3year <- (last_year-(2:0)) %>% as.character()
  last_year <- last_year %>% as.character()

  # define waa catch
  res_sam$input$dat$waa.catch <- waa_catch


  if(res_sam$input$SR=="BHS"){
    res_sam$input$SR <- res_sam$SR <- "HS"
  }

  # define state variables
  aveF  <- colSums(res_sam$caa * res_sam$input$dat$waa.catch * res_sam$faa)/colSums(res_sam$caa * res_sam$input$dat$waa)
  Erate <- colSums(res_sam$caa * res_sam$input$dat$waa.catch)/colSums(res_sam$baa)
  ssb   <- colSums(res_sam$ssb)

  # calculate biological reference points
  Fcurrent <- rowMeans(res_sam$faa[, as.character(year_Fcur)])
  refs <- frasyr::ref.F(res_sam, Fcurrent=Fcurrent, pSPR=perSPR,rps.year=as.numeric(colnames(res_sam$naa)), 
                        M.year=year_biol, waa.year=year_biol, maa.year=year_biol,plot=FALSE)
  if(refs$Fmed[1]==0){
    refs <- frasyr::ref.F(res_sam, Fcurrent=Fcurrent, pSPR=perSPR,rps.year=as.numeric(colnames(res_sam$naa)), 
                          M.year=year_biol, waa.year=year_biol, maa.year=year_biol,plot=FALSE, Fem.init=-2) ## ad hoc fix  
  }
  currentSPR <- refs$currentSPR$perSPR
  refs <- refs$summary
  refs <- refs[colnames(refs)%in%c("Fmed","F0.1",str_c("FpSPR.",perSPR,".SPR"))][3,]
  colnames(refs) <- str_c(colnames(refs),"/Fcur")
    

  # calculate MSY reference points (use function copied from OMutility, OMutilityのときは6歳のF=1と定義していた。ここではどう定義する？)
  biopar <- derive_biopar(res_sam, year_biol)
  if(res_sam$input$SR!="Prop"){
    logb <- res_sam$par_list$rec_logb
    if(is.null(logb)) logb <- res_sam$input$p0.list$rec_logb
    RFmsy   <- Calcu_Fmsy(M=biopar$M, Sel=Fcurrent, w=biopar$waa, g=biopar$maa, method="Baranov", alpha=exp(res_sam$par_list$rec_loga), beta=exp(logb), method_SR=res_sam$SR)
    Bmsy   <- Calcu_Bmsy(RFmsy, M=biopar$M, Sel=Fcurrent, w=biopar$waa, g=biopar$maa, alpha=exp(res_sam$par_list$rec_loga), beta=exp(logb), method_SR=res_sam$SR)
    SBmsy   <- Calcu_SBmsy(RFmsy, M=biopar$M, Sel=Fcurrent, w=biopar$waa, g=biopar$maa, alpha=exp(res_sam$par_list$rec_loga), beta=exp(logb), method_SR=res_sam$SR)

    # calculate additional reference values
    refs_msy <- frasyr::ref.F(res_sam, Fcurrent=Fcurrent*RFmsy, pSPR=perSPR,rps.year=as.numeric(colnames(res_sam$naa)), 
                              M.year=year_biol, waa.year=year_biol, maa.year=year_biol,plot=FALSE)
    FmsySPR <- refs_msy$currentSPR$perSPR
    aa <- calc_steepness(SR=res_sam$SR, rec_pars=make_SRres(res_sam)$pars, M=biopar$M, waa=biopar$waa, maa=biopar$maa,Pope=FALSE,faa=Fcurrent)
    #  expect_equal(aa$Bmsy,Bmsy, tol=0.0001)
    # expect_equal(aa$Fmsy2F, RFmsy, tol=1e-4)
    h <- aa$h
    SB0 <- aa$SB0
  }
  else{
    RFmsy <- Bmsy <- SBmsy <- h <- SB0 <- FmsySPR <- NA
  }
  TBy <- res_sam$baa[,last_year] %>% sum()
  SBy <- res_sam$ssb[,last_year] %>% sum()
    
  res <- tribble(~stat, ~value,
                 str_c("TBy",last_year), TBy,
                 str_c("Sby",last_year), SBy,
                 str_c("Ry", last5year), res_sam$naa[1,last5year],
                 str_c("AFy",last5year), aveF[last5year],
                 str_c("Ey", last5year), Erate[last5year],
                 "currentSPR"  , currentSPR,                 
                 str_c("deple_median_last3"), mean(ssb[last3year])/median(ssb),
                 names(refs) , as.numeric(unlist(refs)),
                 "RFmsy"     , RFmsy,
                 "Bmsy"      , Bmsy,
                 "SBmsy"     , SBmsy,
                 "h"         , h,
                 "SB0"       , SB0,
                 "SBmsy/SB0" , SBmsy/SB0,
                 "FmsySPR"   , FmsySPR,
                 "B/Bmsy"     , TBy/Bmsy/1000,
                 "SB/SBmsy"    , SBy/SBmsy/1000,
                 "SBmsy/SBmax"    , SBmsy*1000/max(ssb)
                 ) %>%
    unnest(cols=c(stat, value))

  return(res)
}

#' @export
#' 
#' @encoding UTF-8
#' 

make_SRres <- function(res_sam, multi_ssb=100, get_rand=FALSE, nsim=10000){
  res_sam2 <- res_sam
  res_sam2$input$Pope <- FALSE
  res_sam2$rec.par <- c(res_sam2$rec.par, sd=res_sam2$sigma.logN[1],rho=0)
  res_SR <- list(pars=as.list(res_sam2$rec.par))
  res_SR$input <- res_sam2$input
  res_SR$input$type <- "L2"
  res_SR$input$SRdata <- get.SRdata(res_sam)
  res_SR$input$SR <- res_sam$input$SR
  ssb <- seq(from=1,to=multi_ssb*max(colSums(res_sam2$ssb)),length=100)
  if(res_SR$input$SR=="BHS"){
    res_SR$pred <- tibble(SSB=ssb/1000,R=1000*frasyr::SRF_BHS(ssb/1000,res_SR$pars$a,res_SR$pars$b,1))
  }
  if(res_SR$input$SR=="BH"){
    res_SR$pred <- tibble(SSB=ssb/1000,R=1000*frasyr::SRF_BH(ssb/1000,res_SR$pars$a,res_SR$pars$b,1))
  }
  res_SR$AICc <- NA

  return(res_SR)
}

#' @export
#' 
#' @encoding UTF-8
#' 

get_randpar <- function(res_sam, nsim=10000, only_fix=TRUE){
  prec1 = res_sam$rep$jointPrecision
  covar1 <- solve(prec1)

# definition of rmvnorm_prec
  rmvnorm_prec <- function(mu, prec, n.sims,seed=1) {
    set.seed(seed)
    z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L_inv = Matrix::Cholesky(prec, super=TRUE)
    return(mu + solve(as(L_inv, 'pMatrix'), solve(t(as.matrix(as(L_inv, 'Matrix'))), z)))
  }

  if(only_fix==TRUE){
    fix.par <- res_sam$opt$par
    coname <- names(res_sam$opt$par)
    fix.prec <- prec1[1:length(res_sam$opt$par),1:length(res_sam$opt$par)]
    fix.covar <- covar1[1:length(res_sam$opt$par),1:length(res_sam$opt$par)]
    parname <- coname[1:length(res_sam$opt$par)]
    
#    rand1 <- rmvnorm_prec(fix.par,fix.prec,n.sims=nsim)
    rand <- MASS::mvrnorm(n=nsim,mu=fix.par,Sigma=fix.covar)  
  }
  else{
    mu_est <- res_sam$obj$env$last.par.best
    rand <- MASS::mvrnorm(n=nsim,mu=mu_est,Sigma=covar1) 
  }
  rand %>% as_tibble(.name_repair="unique") %>% return()
}

#' @export
#' 
#' @encoding UTF-8
#'
#' @examples
#' 

make_samrand <- function(res_sam, nsim=1000){

  rand_par <- get_randpar(res_sam, nsim=nsim, only_fix=FALSE)
  naa_faa <- rand_par %>% select(starts_with("U")) %>%
    unlist() %>% array(dim=c(nsim, 13, 53))
  res_rand <- list()
  for(i in 1:nsim){
    res_rand[[i]] <- res_sam
    res_rand[[i]]$par_list$rec_loga <- rand_par$rec_loga[i]
    res_rand[[i]]$par_list$rec_logb <- rand_par$rec_logb[i]
    if(!is.null(rand_par$rec_loga)) res_rand[[i]]$rec.par["a"] <- rand_par$rec_loga[i] %>% exp()
    if(!is.null(rand_par$rec_logb)) res_rand[[i]]$rec.par["b"] <- rand_par$rec_logb[i] %>% exp()    
    res_rand[[i]]$sigma.logN[1] <- rand_par$logSdLogN[i] %>% exp()
    res_rand[[i]]$naa[] <- naa_faa[i,1:7,] %>% exp()
    res_rand[[i]]$faa[] <- naa_faa[i,c(8:13,13),] %>% exp()
    res_rand[[i]]$ssb <-  res_rand[[i]]$naa * res_rand[[i]]$input$dat$waa * res_rand[[i]]$input$dat$maa
    res_rand[[i]]$baa <-  res_rand[[i]]$naa * res_rand[[i]]$input$dat$waa
    res_rand[[i]]$caa <-  catch_equation(naa=res_rand[[i]]$naa, M=res_rand[[i]]$input$dat$M, faa=res_rand[[i]]$faa, waa=1, Pope=0)
  }
  return(res_rand)
}

#'
#' @export
#' 

# calculate all confidence interval
do_allboot <- function(res_sam_list, nsim=100, year_biol=2020:2022){
    res_pm <- NULL
    res_SR_boot <- res_SR_fit <- res_sam_boot <- list()
  for(i in 1:length(res_sam_list)){
    res_sam_boot[[i]] <- make_samrand(res_sam_list[[i]],nsim=nsim)
    res_pm_boot <- purrr::map_dfr(res_sam_boot[[i]], function(x) get_pm(x,year_biol=year_biol), .id="sim") %>%
        mutate(Model=i)
    res_pm <- bind_rows(res_pm,res_pm_boot)
    res_SR_boot[[i]] <- purrr::map(res_sam_boot[[i]], function(x) make_SRres(x, multi_ssb=1))
    res_SR_fit[[i]] <- purrr::map(res_sam_boot[[i]], function(x){
        SRdata <- get.SRdata(x)
        fit.SR(SRdata, SR="BH")
    })
  }
  list(pm=res_pm, SR=res_SR_boot, SR_fit=res_SR_fit, sam=res_sam_boot)
}
