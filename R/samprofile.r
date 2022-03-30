

#' SAMで推定されたパラメータに対してプロファイル尤度を計算する
#'
#' @param samres sam object
#' @param param_name 変化させるパラメータの名前
#' @param which_param 変化させるパラメータの位置（複数のパラメータが同じ名前を持つときに使用）
#'
#' @export

samprofile <- function(samres,param_name=NULL,which_param = NULL,param_range=NULL,length=50) {

  data = samres$data
  init = samres$par_list
  map = samres$map
  map[["rec_logb"]] <- factor(NA)

  if (is.null(param_range)) {
    param_est = samres$par_list
  }

  if (is.null(param_name) & is.null(which_param)) {
    stop("Either name or linconb must be used")
  }
  if (is.null(which_param)) {
    if (length(init[[param_name]])!=0) stop("The augument of 'param_name' not appropriate. Use 'which_param' if the number of parameters having 'param_name' is more than one")
    init[[param_name]] <- factor(NA)
  } else {
    names(samres$opt$par)[which_param]
  }


}
