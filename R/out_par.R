
#' SAMのFixed effect parametersの表を出力する
#'
#' @export
#'
out_par = function(res,filename=NULL) {
  param_table = data.frame("FE"=names(res$rep$par.fixed),"MLE" = res$rep$par.fixed,"SE"=sqrt(diag(res$rep$cov.fixed))) %>%
    mutate("Gradient" = res$rep$gradient.fixed) %>%
    mutate("Unlinked value" = ifelse(str_detect(FE,"logit"),1/(1+exp(-MLE)),exp(MLE)))

  if(is.null(filename)) filename <- "FE_parameter_table"
  write.csv(param_table,file=paste0(filename,".csv"),row.names = FALSE)
}
