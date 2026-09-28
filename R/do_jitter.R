
#' Do jitter analysis
#'
#' @param samres sam object
#' @param SD SD for random normal distribution to be added to initial values
#' @param nsim default: 100
#'
#' @export
#'
do_jitter = function(samres,SD=0.1,nsim=100,seed=1) {
  res <- samres
  init <- init0 <- res$init
  map <- res$map
  # map$logB
  parname = unique(names(res$opt$par))
  tmbdata <- res$data
  random <- res$obj$env$random
  cpp.file.name <- res$input$cpp.file.name
  set.seed(seed)
  cat(paste0("0: ",round(res$opt$objective,3),"\n"))
  resdat = data.frame(ID=0,obj_value=res$opt$objective)
  reslist = list()
  for(k in 1:nsim) {
    for(j in 1:length(parname)) {
      if (is.null(map[[parname[j]]])) {
        init[[parname[j]]] <- rnorm(length(init0[[parname[j]]]),init0[[parname[j]]],SD)
      } else {
        loc = which(!is.na(map[[parname[j]]]),TRUE)
        init[[parname[j]]][loc] <- rnorm(length(loc),init0[[parname[j]]][loc],SD)
      }
    }
    # input <- res$input

    tmp = try(est_fixed(tmbdata,par_init=init,map=map,random=random,cpp.file.name = cpp.file.name,silent=TRUE))

    if(class(tmp)[1] != "try-error") {
      if (res$opt$objective-tmp$opt$objective > 0) {
        init0 <- init
        res <- tmp
      }
    }else{
      tmp$opt$objective <- NA
    }
    resdat = bind_rows(resdat,data.frame(ID=k,obj_value=tmp$opt$objective))
    reslist[[k]] <- tmp
    cat(paste0(k,": ", round(tmp$opt$objective,3),"\n"))
  }
  return = list(resdat=resdat,reslist=reslist)
  return( return )
}
