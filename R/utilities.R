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
                   res_future,
                   waa_catch,                     
                   last_year=2022,
                   waa_biom=res_sam$input$dat$waa,
                   year_biol=2020:2022,
                   year_Fcur=2020:2022,
                   perSPR=c(30,40,50)
                   ){

  # define year range
  last5year <- (last_year-(4:0)) %>% as.character()
  last3year <- (last_year-(2:0)) %>% as.character()
  last_year <- last_year %>% as.character()

  # define waa catch
  res_sam$input$dat$waa.catch <- waa_catch 

  # define state variables
  aveF  <- colSums(res_sam$caa * res_sam$input$dat$waa.catch * res_sam$faa)/colSums(res_sam$caa * res_sam$input$dat$waa)
  Erate <- colSums(res_sam$caa * res_sam$input$dat$waa.catch)/colSums(res_sam$baa)
  ssb   <- colSums(res_sam$ssb)

  # calculate biological reference points
  Fcurrent <- rowMeans(res_sam$faa[, as.character(year_Fcur)])
  refs <- frasyr::ref.F(res_sam, Fcurrent=Fcurrent, pSPR=perSPR,
                        M.year=year_biol, waa.year=year_biol, maa.year=year_biol,plot=FALSE)$summary
  refs <- refs[colnames(refs)%in%c("F0.1",str_c("FpSPR.",perSPR,".SPR"))][3,]

  # calculate MSY reference points (use function copied from OMutility, OMutilityのときは6歳のF=1と定義していた。ここではどう定義する？)
  biopar <- derive_biopar(res_sam, year_biol)
  RFmsy   <- Calcu_Fmsy(M=biopar$M, Sel=biopar$faa, w=biopar$waa, g=biopar$maa, method="Baranov", alpha=exp(res_sam$par_list$rec_loga), beta=exp(res_sam$par_list$rec_logb), method_SR=res_sam$SR)
  Bmsy   <- Calcu_Bmsy(RFmsy, M=biopar$M, Sel=biopar$faa, w=biopar$waa, g=biopar$maa, alpha=exp(res_sam$par_list$rec_loga), beta=exp(res_sam$par_list$rec_logb), method_SR=res_sam$SR)
  SBmsy   <- Calcu_SBmsy(RFmsy, M=biopar$M, Sel=biopar$faa, w=biopar$waa, g=biopar$maa, alpha=exp(res_sam$par_list$rec_loga), beta=exp(res_sam$par_list$rec_logb), method_SR=res_sam$SR)
  TBy <- res_sam$baa[,last_year] %>% sum()
  SBy <- res_sam$ssb[,last_year] %>% sum()
    
  res <- tribble(~stat, ~value,
                 str_c("TBy",last_year), TBy,
                 str_c("Sby",last_year), SBy,
                 str_c("Ry", last5year), res_sam$naa[1,last5year],
                 str_c("AFy",last5year), aveF[last5year],
                 str_c("Ey", last5year), Erate[last5year],
                 str_c("deple_median_last3"), mean(ssb[last3year])/median(ssb),
                 names(refs) , as.numeric(unlist(refs)),
                 "RFmsy"     , RFmsy,
                 "Bmsy"      , Bmsy,
                 "SBmsy"     , SBmsy,
                 "RBmsy"     , TBy/Bmsy/1000,
                 "RSBmsy"    , SBy/SBmsy/1000
                 ) %>%
    unnest(cols=c(stat, value))

  return(res)
}
