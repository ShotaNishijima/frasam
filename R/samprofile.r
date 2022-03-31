

#' SAMで推定されたパラメータに対してプロファイル尤度を計算する
#'
#' @param samres sam object
#' @param param_name 変化させるパラメータの名前
#' @param which_param 変化させるパラメータの位置（複数のパラメータが同じ名前を持つときに使用）
#'
#' @export

samprofile <- function(samres,param_name,which_param=1,param_range=NULL,length=50) {

  data = samres$data
  init = samres$par_list
  map = samres$map

  if (is.null(param_range)) {
    param_est = samres$par_list[[param_name]][which_param]
    param_range = param_est + c(-2,2)
  }

  prof_param = seq(param_range[1],param_range[2],length=length)

  map[[param_name]] <- 0:(length(init[[param_name]])-1)
  map[[param_name]][which_param] <- NA
  map[[param_name]] <- factor(map[[param_name]])

  random = unique(names(samres$rep$par.random))

  obj <- opt <- list()
  obj_tbl = tibble(par=samres$par_list[[param_name]][which_param],
                   obj_value = samres$opt$objective,
                   likelihood="maximum")

  par_tbl = tibble(name = names(samres$opt$par),value = samres$opt$par,
                   likelihood="maximum")

  for (i in 1:length) {
    init[[param_name]][which_param] <- prof_param[i]
    obj_tmp = TMB::MakeADFun(data,init,map=map,random=random,DLL=samres$input$cpp.file.name,silent=TRUE)
    opt_tmp = nlminb(obj_tmp$par, obj_tmp$fn, obj_tmp$gr)
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
