
#' Update SAM result by new parameters of fixed and random effects
#'
#' @param samres SAM object
#' @param new_par new parameter set
#'
#' @export
#'
update_sam <- function(samres,new_par) {
  obj_orig <- obj_update <- samres$obj
  if (length(obj_update$env$last.par) != length(new_par)) {
    stop("'new_par' has a different length from 'samres$obj$env$last.par'")
  } else {
    obj_update$env$last.par <- obj_update$env$last.par.best <- new_par
  }
  fixed_name <- unique(names(obj_orig$par))
  fixed_sim <- obj_orig$par
  if (is.null(names(new_par))) stop("'names(new_par)' is required")
  for (i in 1:length(fixed_name)) {
    fixed_sim[names(obj_orig$par) == fixed_name[i]] <- new_par[names(new_par)==fixed_name[i]]
  }
  # cbind(fixed_sim, obj_orig$par)
  obj_update$par <- fixed_sim

  input = samres$input
  input$obj_overwrite <- obj_update
  res = safe_do_call(sam,input)

  return( res )
}
