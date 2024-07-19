
#' Do leave-one-out index analysis
#'
#' Indexを一つずつ除いて解析するための関数
#'
#' @inheritParams retro_sam
#' @export
do_loo_index = function(samres) {
  res <- samres
  input = res$input
  input$p0.list <- NULL
  nindex = nrow(res$input$dat$index)

  RES = lapply(1:nindex,function(j) {
    input$dat$index <- res$input$dat$index[-j,]
    input$a.init <- as.numeric(res$rec.par["a"])
    input$b.init <- as.numeric(res$rec.par["b"])
    input$q.init <- res$q[-j]
    input$abund <- res$input$abund[-j]
    input$min.age <- res$input$min.age[-j]
    input$max.age <- res$input$max.age[-j]
    if(!is.null(res$input$b.fix)) input$b.fix <- res$input$b.fix[-j]
    if(!is.null(res$input$index.key)) input$index.key <- res$input$index.key[-j] - min(res$input$index.key[-j])
    if(!is.null(res$input$index.b.key)) input$index.b.key <- res$input$index.b.key[-j] - min(res$input$index.b.key[-j])
    if(!is.null(res$input$catch_prop)) input$catch_prop <- res$input$catch_prop[,,-j]

    input$no_est <- TRUE
    tmp = do.call(sam,input)

    p0_list <- tmp$obj$env$parList()
    p0_list$U <- res$obj$env$parList()[["U"]]
    p0_list$logSdLogFsta <- res$obj$env$parList()[["logSdLogFsta"]]
    p0_list$logSdLogN <- res$par_list$logSdLogN

    input$no_est <- FALSE
    input$p0.list <- p0_list

    res.c = do.call(sam,input)
    return( res.c )
  })
  return( RES )
}
