#'
#' SAMの結果を読み込んで主要なパラメータの予測値とSD、信頼区間などを出力する
#' 
#' @export
#' 

make_assess_result <- function(res=res, CI=0.90){
  vcvdata <- tibble(stat0=names(res$rep$value), Value=res$rep$value, SD=res$rep$sd, CV=res$rep$sd/res$rep$value)
  vcvdata <- vcvdata %>% mutate(Cz = exp(qnorm(CI+(1-CI)/2)*sqrt(log(1+CV^2)))) %>% mutate(lower = Value/Cz, upper = Value*Cz)
  vcvdata$Obs <- NA
  vcvdata$Year <- NA
  vcvdata$Age <- NA
  vcvdata[vcvdata$stat0=="exp_logN",]$Age <-  rep(as.numeric(rownames(res$naa)),ncol(res$naa))
  vcvdata[vcvdata$stat0=="exp_logN",]$Year <- sort(rep(as.numeric(colnames(res$naa)), length.out=ncol(res$naa)*nrow(res$naa)))
  vcvdata[vcvdata$stat0=="exp_logF",]$Age <-  rep(as.numeric(rev(rev(rownames(res$naa))[-1])),ncol(res$naa))
  vcvdata[vcvdata$stat0=="exp_logF",]$Year <- sort(rep(as.numeric(colnames(res$naa)), length.out=ncol(res$naa)*(nrow(res$naa)-1)))
  
  vcvdata[vcvdata$stat0=="ssb",]$Year <- (as.numeric(colnames(res$naa)))
  vcvdata[vcvdata$stat0=="B_total",]$Year <- (as.numeric(colnames(res$naa)))
  vcvdata[vcvdata$stat0=="F_mean",]$Year <- (as.numeric(colnames(res$naa)))
  if(is.null(res$input$dat$waa.catch)){
    res$input$dat$waa.catch <- res$input$dat$waa
  }
  wcaa <- res$input$dat$caa * res$input$dat$waa.catch
  vcvdata[vcvdata$stat0=="Catch_biomass",]$Obs <- (as.numeric(colSums(wcaa)))
  vcvdata[vcvdata$stat0=="Catch_biomass",]$Year <- (as.numeric(colnames(res$naa)))
  vcvdata[vcvdata$stat0=="Exploitation_rate",]$Year <- (as.numeric(colnames(res$naa)))
  
  vcvdata[vcvdata$stat0=="stockMeanWeight_true",]$Age <-  rep(as.numeric(rownames(res$naa)),ncol(res$naa))
  vcvdata[vcvdata$stat0=="stockMeanWeight_true",]$Year <- sort(rep(as.numeric(colnames(res$naa)), length.out=ncol(res$naa)*nrow(res$naa)))
  
  vcvdata <- vcvdata %>% mutate(stat = case_when(stat0=="ssb" ~ "SSB",
                                                 stat0=="F_mean" ~ "F",
                                                 stat0=="B_total" ~ "Biomass",
                                                 stat0=="Catch_biomass" ~ "Catch",
                                                 stat0=="exp_logF" ~ "faa",
                                                 stat0=="exp_logN" ~ "naa",
                                                 stat0=="Exploitation_rate" ~ "FishingRatio",
                                                 stat0=="stockMeanWeight_true" ~ "waa",
                                                 TRUE~stat0))
  vcvdata$Model <- "SAM"
  rec <- vcvdata[vcvdata$stat0=="exp_logN" & vcvdata$Age==min(vcvdata$Age,na.rm=T),] 
  rec$stat <- rec$stat0 <- "Recruitment"
  vcvdata <- rbind(vcvdata, rec)
  vcvdata
}
