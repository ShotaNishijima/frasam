#' Check convergence diagnostics for a SAM fit
#'
#' \code{check_fit_sam} summarizes basic convergence diagnostics for a
#' \code{sam()} result object. The function does not change the fitted object.
#'
#' @param res A \code{sam} result object.
#' @param gradient_tol Maximum allowed absolute fixed-effect gradient.
#' @param se_max Maximum allowed fixed-effect standard error. Set \code{Inf} to
#'   skip this check.
#' @param sigma_range Allowed range for reported standard deviations.
#' @param rho_range Allowed range for reported rho values.
#' @param par_abs_max Maximum allowed absolute fixed-effect estimate on the
#'   internal scale. Set \code{Inf} to skip this check.
#' @param boundary_tol Tolerance for detecting parameters close to finite
#'   lower or upper bounds.
#' @param verbose If \code{TRUE}, print the diagnostic table.
#'
#' @return A list with \code{ok}, \code{checks}, \code{fixed}, \code{sigma},
#'   and \code{boundary} elements.
#' @export
check_fit_sam <- function(res,
                          gradient_tol = 1e-2,
                          se_max = Inf,
                          sigma_range = c(1e-4, 10),
                          rho_range = c(1e-4, 1 - 1e-4),
                          par_abs_max = Inf,
                          boundary_tol = 1e-4,
                          verbose = TRUE) {
  add_check <- function(checks, check, ok, value = NA, threshold = NA, message = "") {
    rbind(
      checks,
      data.frame(
        check = check,
        ok = as.logical(ok),
        value = as.character(value),
        threshold = as.character(threshold),
        message = message,
        stringsAsFactors = FALSE
      )
    )
  }

  checks <- data.frame(
    check = character(),
    ok = logical(),
    value = character(),
    threshold = character(),
    message = character(),
    stringsAsFactors = FALSE
  )

  if (!inherits(res, "sam")) {
    stop("'res' must be a sam result object.", call. = FALSE)
  }

  opt <- res$opt
  rep <- res$rep

  convergence <- if (!is.null(opt$convergence)) opt$convergence else NA_integer_
  checks <- add_check(
    checks,
    "optimizer convergence",
    isTRUE(convergence == 0),
    convergence,
    "0",
    if (isTRUE(convergence == 0)) "nlminb convergence code is 0." else "nlminb convergence code is not 0."
  )

  pdHess <- if (!is.null(rep$pdHess)) rep$pdHess else NA
  checks <- add_check(
    checks,
    "positive definite Hessian",
    isTRUE(pdHess),
    pdHess,
    "TRUE",
    if (isTRUE(pdHess)) "sdreport reports pdHess = TRUE." else "sdreport does not report pdHess = TRUE."
  )

  gradient <- if (!is.null(rep$gradient.fixed)) rep$gradient.fixed else NA_real_
  max_gradient <- suppressWarnings(max(abs(gradient), na.rm = TRUE))
  if (!is.finite(max_gradient)) max_gradient <- NA_real_
  checks <- add_check(
    checks,
    "maximum fixed-effect gradient",
    is.finite(max_gradient) && max_gradient <= gradient_tol,
    signif(max_gradient, 4),
    gradient_tol,
    "Maximum absolute fixed-effect gradient."
  )

  fixed_est <- if (!is.null(rep$par.fixed)) rep$par.fixed else numeric()
  fixed_se <- if (!is.null(rep$cov.fixed)) sqrt(diag(rep$cov.fixed)) else rep(NA_real_, length(fixed_est))
  if (length(fixed_se) != length(fixed_est)) fixed_se <- rep(NA_real_, length(fixed_est))
  if (is.null(names(fixed_se))) names(fixed_se) <- names(fixed_est)
  fixed_names <- names(fixed_est)
  if (is.null(fixed_names)) fixed_names <- paste0("par", seq_along(fixed_est))

  fixed <- data.frame(
    name = fixed_names,
    estimate = as.numeric(fixed_est),
    se = as.numeric(fixed_se),
    gradient = as.numeric(gradient[seq_along(fixed_est)]),
    stringsAsFactors = FALSE
  )

  finite_fixed <- is.finite(fixed$estimate) & is.finite(fixed$se)
  checks <- add_check(
    checks,
    "finite fixed effects and SE",
    all(finite_fixed),
    paste0(sum(finite_fixed), "/", nrow(fixed)),
    "all",
    "Fixed-effect estimates and standard errors should be finite."
  )

  max_se <- suppressWarnings(max(fixed$se, na.rm = TRUE))
  if (!is.finite(max_se)) max_se <- NA_real_
  checks <- add_check(
    checks,
    "maximum fixed-effect SE",
    is.infinite(se_max) || (is.finite(max_se) && max_se <= se_max),
    signif(max_se, 4),
    se_max,
    "Large standard errors can indicate weak identification."
  )

  max_abs_par <- suppressWarnings(max(abs(fixed$estimate), na.rm = TRUE))
  if (!is.finite(max_abs_par)) max_abs_par <- NA_real_
  checks <- add_check(
    checks,
    "maximum absolute fixed effect",
    is.infinite(par_abs_max) || (is.finite(max_abs_par) && max_abs_par <= par_abs_max),
    signif(max_abs_par, 4),
    par_abs_max,
    "Extremely large internal-scale estimates can indicate weak identification."
  )

  sigma_names <- c("sigma", "sigma.logC", "sigma.logFsta", "sigma.logN")
  sigma <- do.call(
    rbind,
    lapply(sigma_names, function(nm) {
      x <- res[[nm]]
      if (is.null(x)) return(NULL)
      data.frame(
        type = nm,
        index = seq_along(x),
        value = as.numeric(x),
        stringsAsFactors = FALSE
      )
    })
  )
  if (is.null(sigma)) {
    sigma <- data.frame(type = character(), index = integer(), value = numeric(), stringsAsFactors = FALSE)
  }
  sigma$ok <- is.finite(sigma$value) & sigma$value >= sigma_range[1] & sigma$value <= sigma_range[2]
  sigma$problem <- ifelse(
    is.na(sigma$value) | !is.finite(sigma$value),
    "non-finite",
    ifelse(sigma$value < sigma_range[1], "too small", ifelse(sigma$value > sigma_range[2], "too large", ""))
  )
  sigma_ok <- nrow(sigma) == 0 || all(sigma$ok)
  checks <- add_check(
    checks,
    "reported sigma range",
    sigma_ok,
    if (nrow(sigma) == 0) NA else paste(signif(range(sigma$value, na.rm = TRUE), 4), collapse = " to "),
    paste(sigma_range, collapse = " to "),
    "Reported standard deviations should be finite and within the diagnostic range."
  )

  rho_values <- if (!is.null(res$rho)) res$rho else numeric()
  rho_values <- rho_values[!is.na(rho_values)]
  rho_fixed_at_boundary <- !is.null(res$input$rho.mode) && res$input$rho.mode %in% c(0, 1)
  rho_ok <- isTRUE(rho_fixed_at_boundary) ||
    length(rho_values) == 0 ||
    all(is.finite(rho_values) & rho_values >= rho_range[1] & rho_values <= rho_range[2])
  checks <- add_check(
    checks,
    "reported rho range",
    rho_ok,
    if (length(rho_values) == 0) NA else paste(signif(range(rho_values, na.rm = TRUE), 4), collapse = " to "),
    paste(rho_range, collapse = " to "),
    if (isTRUE(rho_fixed_at_boundary)) {
      "Rho is fixed by rho.mode, so the boundary range check is skipped."
    } else {
      "Rho values very close to 0 or 1 can indicate boundary behavior."
    }
  )

  boundary <- data.frame(
    name = character(),
    estimate = numeric(),
    lower = numeric(),
    upper = numeric(),
    near_lower = logical(),
    near_upper = logical(),
    stringsAsFactors = FALSE
  )

  lower <- res$input$lower
  upper <- res$input$upper
  opt_par <- if (!is.null(opt$par)) opt$par else numeric()
  if (!is.null(lower) || !is.null(upper)) {
    if (is.null(lower)) lower <- rep(-Inf, length(opt_par))
    if (is.null(upper)) upper <- rep(Inf, length(opt_par))
    if (length(lower) == length(opt_par) && length(upper) == length(opt_par)) {
      opt_names <- names(opt_par)
      if (is.null(opt_names)) opt_names <- paste0("par", seq_along(opt_par))
      boundary <- data.frame(
        name = opt_names,
        estimate = as.numeric(opt_par),
        lower = as.numeric(lower),
        upper = as.numeric(upper),
        near_lower = is.finite(lower) & abs(opt_par - lower) <= boundary_tol,
        near_upper = is.finite(upper) & abs(opt_par - upper) <= boundary_tol,
        stringsAsFactors = FALSE
      )
    }
  }
  n_boundary <- sum(boundary$near_lower | boundary$near_upper, na.rm = TRUE)
  checks <- add_check(
    checks,
    "parameters near finite bounds",
    n_boundary == 0,
    n_boundary,
    "0",
    "Only checked when finite lower or upper bounds were supplied."
  )

  ok <- all(checks$ok)
  out <- list(
    ok = ok,
    checks = checks,
    fixed = fixed,
    sigma = sigma,
    boundary = boundary
  )
  class(out) <- "check_fit_sam"

  if (isTRUE(verbose)) {
    print(checks, row.names = FALSE)
    if (!isTRUE(ok)) warning("Some convergence diagnostics failed.", call. = FALSE)
  }

  invisible(out)
}
