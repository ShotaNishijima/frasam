

#' SAMで推定されたパラメータに対してプロファイル尤度を計算する
#'
#' @param samres sam object
#' @param param_name 変化させるパラメータの名前
#' @param which_param 変化させるパラメータの位置（複数のパラメータが同じ名前を持つときに使用）
#'
#' @encoding UTF-8
#'
#' @export

samprofile <- function(samres,param_name,which_param=1,param_range=NULL,length=50) {

  data = samres$data
  init = samres$par_list
  # init = samres$opt$par
  map = samres$map

  if (is.null(param_range)) {
    param_est = samres$par_list[[param_name]][which_param]
    param_range = param_est + c(-2,2)
  }

  prof_param = seq(param_range[1],param_range[2],length=length)

  if(is.null(map[[param_name]])) {
    map[[param_name]] <- 0:(length(samres$par_list[[param_name]])-1)
  }
  map[[param_name]][which_param] <- NA
  map[[param_name]] <- factor(map[[param_name]])

  random = unique(names(samres$rep$par.random))

  obj <- opt <- list()
  obj_tbl = tibble(par=samres$par_list[[param_name]][which_param],
                   obj_value = samres$opt$objective,
                   likelihood="maximum")

  par_tbl = tibble(name = names(samres$opt$par),value = samres$opt$par,
                   likelihood="maximum")

  # obj = samres$obj
  # tmp <- try(obj$fn(opj$par))
  # if(class(tmp)[1]=="try-error") {
  # obj = TMB::MakeADFun(data,init,map=map,random=random,DLL=samres$input$cpp.file.name,silent=TRUE)
  # }
  # obj$fn(opj$par)
  nlminb.control = list(eval.max = 1e4,
                        iter.max = 1e4,
                        trace = 0)
  for (i in 1:length) {
    init[[param_name]][which_param] <- prof_param[i]
    obj_tmp = TMB::MakeADFun(data,init,map=map,random=random,DLL=samres$input$cpp.file.name,silent=TRUE)
    obj_tmp$fn(opt_tmp$par)
    opt_tmp = nlminb(obj_tmp$par, obj_tmp$fn, obj_tmp$gr,control=nlminb.control)
    # opt_tmp = nlminb(obj_tmp$par, obj_tmp$fn, obj_tmp$gr,control=nlminb.control)
    message(paste0("par: ",prof_param[i],"   objective: ", opt_tmp$objective))
    obj[[i]] <- obj_tmp
    opt[[i]] <- opt_tmp

    obj_tbl = bind_rows(obj_tbl,
                        tibble(par=prof_param[i],
                               obj_value = opt_tmp$objective,
                               likelihood="profile"))

    par_tbl = par_tbl %>% bind_rows(tibble(name = param_name,value = prof_param[i],
           likelihood="profile")) %>%
      bind_rows(tibble(name = names(opt_tmp$par),value = opt_tmp$par,
                           likelihood="profile")
      )
  }

  return (list(obj_tbl=obj_tbl,par_tbl=par_tbl,obj_list=obj,opt_list=opt))
}


#' SAMの固定効果だけを推定する関数
#'
#' @param tmbdata tmbdata
#' @param par_init パラメータの初期値
#' @param map 推定しない（初期値で固定する）パラメータ
#' @param random ランダム効果で推定するパラメータ
#' @param cpp.file.name cppファイルの名前（初期設定："sam"）
#'
#' @encoding UTF-8
#'
#' @export

est_fixed <- function(tmbdata,par_init,map,random="U",cpp.file.name="sam2",silent=TRUE,...) {

  obj_tmp = TMB::MakeADFun(tmbdata,par_init,map=map,random=random,DLL=cpp.file.name,silent=silent,...)
  obj_tmp$fn(opj_tmp$par)
  nlminb.control = list(eval.max = 1e4,
                        iter.max = 1e4,
                        trace = 0)
  opt_tmp = nlminb(obj_tmp$par, obj_tmp$fn, obj_tmp$gr,control=nlminb.control)

  return(list(obj=obj_tmp,opt=opt_tmp))
}

#' TMB::MakeADFunに必要な引数から固定効果とランダム効果を推定する関数
#'
#' @inheritParams est_fixed
#' @param bias_correct 平均値のバイアス補正をするかどうか
#'
#' @encoding UTF-8
#'
#' @export

est_mixed <- function(tmbdata,par_init,map,random="U",cpp.file.name="sam2",silent=TRUE,bias_correct=TRUE) {

  tmp = est_fixed(tmbdata,par_init,map=map,random=random,cpp.file.name=cpp.file.name,silent=silent)
  rep = TMB::sdreport(tmp$obj,bias.correct=bias_correct)

  return(list(obj=tmp$obj,opt=tmp$opt,rep=rep))
}
