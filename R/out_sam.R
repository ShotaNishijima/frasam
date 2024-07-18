
#' output sam object results
#' @param res sam object
#'
#' @export

out_sam <- function(res,
                   filename="sam",
                   csvname=NULL) {
  # old.par <- par()
  # exit.func <- function(){
  #   dev.off()
  #   options(warn=0)
  # }
  # on.exit(exit.func())
  #
  if(!is.null(filename)){
    csvname <- paste(filename,".csv",sep="")
    # pdfname <- paste(filename,".pdf",sep="")
  }

  # pdf(pdfname)
  # par(mfrow=c(3,2),mar=c(3,3,2,1))
  options(warn=-1)

  write.table2 <- function(x,title.tmp="",is.plot=FALSE,...){
    if(is.plot){
      if(!is.null(dim(x))){
        matplot(colnames(x),t(x),type="b",ylim=c(0,max(x,na.rm=T)),pch=substr(rownames(x),1,1))
      }
      else{
        barplot(x)
      }
      title(title.tmp)
    }
    if(!is.null(dim(x))){
      tmp <- matrix("",nrow(x)+1,ncol(x)+1)
      tmp[-1,-1] <- as.character(unlist(x))
      tmp[-1,1] <- rownames(x)
      tmp[1,-1] <- colnames(x)
    }
    else{
      tmp <- x
    }
    write.table(tmp,append=T,sep=",",quote=FALSE,file=csvname,col.names=F,row.names=F,...)
  }

  pd <- packageDescription("frasam")
  if(is.null(pd$GithubSHA1)) pd$GithubSHA1 <- "local" # fraysrをload_allした場合コミット番号が記録されないため
  write(paste("# frasam(@",pd$GithubSHA1,") outputs at ",date()," & ",getwd(), sep=""),file=csvname)

  write("# SAM results",file=csvname, append=T)

  write("\n# observed catch at age", file=csvname,append=T)
  write.table2(t(res$input$dat$caa),title.tmp="Observed catch at age")

  write("\n# maturity at age", file=csvname,append=T)
  write.table2(t(res$input$dat$maa),title.tmp="Maturity at age")

  write("\n# weight at age for biomass calculation", file=csvname,append=T)
  write.table2(t(res$input$dat$waa),title.tmp="Weight at age (for biomass)")

  if(!is.null(res$input$dat$waa.catch)){
    write("\n# weight at age for catch calculation", file=csvname,append=T)
    write.table2(t(res$input$dat$waa.catch),title.tmp="Weight at age (for catch)")
  }

  write("\n# M at age",file=csvname,append=T)
  write.table2(t(res$input$dat$M),title.tmp="M at age")

  write("\n# fishing mortality at age",file=csvname,append=T)
  write.table2(t(res$faa),title.tmp="F at age")

  write("\n# numbers at age",file=csvname,append=T)
  write.table2(t(res$naa),title.tmp="Numbers at age")

  write("\n# predicted catch at age",file=csvname,append=T)
  write.table2(t(res$caa),title.tmp="Numbers at age")

  sumrep <- summary(res$rep)
  tmp <- sumrep[rownames(sumrep)=="B_total",]
  rownames(tmp) <- colnames(res$naa)
  write("\n# total stock biomass",file=csvname,append=T)
  write.table2(tmp,title.tmp="Total biomass",is.plot=FALSE)

  tmp <- sumrep[rownames(sumrep)=="ssb",]
  rownames(tmp) <- colnames(res$naa)
  write("\n# spawning stock biomass",file=csvname,append=T)
  write.table2(tmp,title.tmp="Spawning stock biomass",is.plot=FALSE)

  # rownames(sumrep) %>% unique
  tmp <- sumrep[rownames(sumrep)=="F_mean",]
  rownames(tmp) <- colnames(res$naa)
  write("\n# mean F",file=csvname,append=T)
  write.table2(tmp,title.tmp="Mean F",is.plot=FALSE)

  tmp <- sumrep[rownames(sumrep)=="Exploitation_rate",]
  rownames(tmp) <- colnames(res$naa)
  write("\n# exploitation rate",file=csvname,append=T)
  write.table2(tmp,title.tmp="Exploitation rate",is.plot=FALSE)

  tmp <- sumrep[rownames(sumrep)=="Catch_biomass",]
  rownames(tmp) <- colnames(res$naa)
  write("\n# catch biomass",file=csvname,append=T)
  write.table2(tmp,title.tmp="Catch biomass",is.plot=FALSE)

  write("\n# YPR & SPR history ",file=csvname,append=T)
  get.SPR(res)$ysdata %>% rownames_to_column(var="year") %>%
    as_tibble() %>%
    dplyr::select(-"F/Ftarget") %>%
    write_csv(path=csvname,append=T, col_names=TRUE)

  write("\n# stock-recruitment relationship",file=csvname,append=T)
  write.table2(cbind("type"=res$input$SR,data.frame(t(res$rec.par)),"sigmaR"=res$sigma.logN[1]),
               title.tmp="Stock-recuitment relationship",is.plot=FALSE)

  tmp <- sumrep[rownames(sumrep)=="exp_logN",]
  tmp <- cbind(expand.grid(age=rownames(res$naa),year=colnames(res$naa)),tmp) %>%
    dplyr::arrange(age,year)
  rownames(tmp) <- NULL
  write("\n# number at age",file=csvname,append=T)
  write.table2(tmp,title.tmp="Number at age",is.plot=FALSE)

  tmp <- sumrep[rownames(sumrep)=="exp_logF",]
  tmp <- cbind(expand.grid(age=rownames(res$faa[-nrow(res$faa),]),year=colnames(res$naa)),tmp) %>%
    dplyr::arrange(age,year)
  rownames(tmp) <- NULL
  write("\n# F at age",file=csvname,append=T)
  write.table2(tmp,title.tmp="F at age",is.plot=FALSE)
  # dev.off()
}
