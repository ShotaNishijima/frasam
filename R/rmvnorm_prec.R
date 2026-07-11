#' Simulate multivariate normal variables given a mean vector and precision matrix
#'
#' @param mu vector of parameter means
#' @param prec joint precision matrix
#' @param n.sims number of draws
#' @param seed seed number
#'
#' @return length(mu) by n.sims matrix of parameter draws
#'
#' @import Matrix
#' @export
rmvnorm_prec <- function(mu, prec, n.sims, seed=123) {
  set.seed(seed)
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L_inv = Matrix::Cholesky(prec, super=TRUE)
  res = mu + solve(as(L_inv, 'pMatrix'), solve(t(as.matrix(as(L_inv, 'Matrix'))), z))
  rownames(res) <- names(mu)
  return( res )
}
