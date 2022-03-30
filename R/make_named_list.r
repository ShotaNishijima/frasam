#' Making a combined list with original list name(s)
#' @param ... list of \code{class("list")}
#' @encoding UTF-8
#' @export
make_named_list <- function(...){
  object_name <- as.character(eval(substitute(alist(...))))
  x <- list(...)
  vnames <- names(x)
  novnames <- !nzchar(vnames)
  if (length(novnames) > 0L){
    vnames[novnames] <- object_name[novnames]
  } else {
    vnames <- object_name
  }
  setNames(x, vnames)
}
