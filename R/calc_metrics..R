#' Compute common regression metrics (R2, RMSE, MAE, RMSRE, MARE, MedBias, MedRelBias)
#'
#' @title Calculate regression/forecast accuracy and bias metrics
#' @description
#' Calculates a set of accuracy (error) and bias metrics between observed values
#' and predictions:
#' \itemize{
#'   \item \strong{R2\_model}: \(1 - \frac{\sum (y - \hat{y})^2}{\sum (y - \bar{y})^2}\)
#'   \item \strong{RMSE}: \(\sqrt{\text{mean}((\hat{y}-y)^2)}\)
#'   \item \strong{MAE}: \(\text{mean}(|\hat{y}-y|)\)
#'   \item \strong{RMSRE}: \(\sqrt{\text{mean}(((\hat{y}-y)/y)^2)}\) (optionally in \%)
#'   \item \strong{MARE}: \(\text{mean}(|\hat{y}-y|/y)\) (optionally in \%)
#'   \item \strong{MedBias}: \(\text{median}(\hat{y}-y)\)
#'   \item \strong{MedRelBias}: \(\text{median}((\hat{y}-y)/y)\) (optionally in \%)
#' }
#'
#' Relative metrics (RMSRE/MARE/MedRelBias) use the observed values \code{y} in the denominator.
#' If zeros occur in \code{y}, control the behavior via \code{rel_zero}:
#' \itemize{
#'   \item \code{"omit"}: drop pairs where \(y=0\) (default for stability)
#'   \item \code{"epsilon"}: replace \(y\) by \(y + \epsilon\) with small \code{epsilon}
#'   \item \code{"none"}: no guard; may yield \code{Inf}/\code{NaN}
#' }
#'
#' @param y Numeric vector of observed values.
#' @param y_pred Numeric vector of predicted values (same length as \code{y}).
#' @param na_rm Logical; if \code{TRUE} (default), remove \code{NA}/\code{NaN}/\code{Inf} after
#'   forming needed quantities for each metric.
#' @param percent Logical; if \code{TRUE} (default), express relative metrics (RMSRE, MARE, MedRelBias)
#'   in percentage (i.e., multiplied by 100). If \code{FALSE}, leave as proportions.
#' @param rel_zero Character; handling of zeros in \code{y} for relative metrics.
#'   One of \code{"omit"}, \code{"epsilon"}, \code{"none"}. See Details.
#' @param epsilon Numeric; small constant added to denominators when \code{rel_zero = "epsilon"}.
#'   Default \code{1e-8}.
#'
#' @return
#' A one-row \code{data.frame} with columns:
#' \code{R2_model, RMSE, MAE, RMSRE, MARE, MedBias, MedRelBias}.
#'
#' @details
#' \itemize{
#' \item \strong{R2\_model} uses the usual in-sample definition with the sample mean of \code{y}.
#'   If \code{var(y) == 0}, \code{R2_model} is returned as \code{NA}.
#' \item Bias-type metrics retain the sign (systematic over-/under-prediction),
#'   while error-type metrics use absolute value or squares (magnitude only).
#' }
#'
#' @examples
#' y      <- c(10, 12, 9, 15, 0, 8)
#' y_pred <- c(11, 11, 8, 14, 0.2, 9)
#'
#' # Default: percent = TRUE, rel_zero = "omit"
#' calc_metrics(y, y_pred)
#'
#' # Use proportions (not %) and guard zeros with epsilon
#' calc_metrics(y, y_pred, percent = FALSE, rel_zero = "epsilon", epsilon = 1e-6)
#'
#' # No guarding (may produce Inf/NaN if y contains zeros)
#' calc_metrics(y, y_pred, rel_zero = "none")
#'
#' @export
calc_metrics <- function(y,
                         y_pred,
                         na_rm   = TRUE,
                         percent = FALSE,
                         rel_zero = c("omit", "epsilon", "none"),
                         epsilon = 1e-8) {
  rel_zero <- match.arg(rel_zero)

  # basic checks
  if (!is.numeric(y) || !is.numeric(y_pred)) {
    stop("`y` and `y_pred` must be numeric vectors.")
  }
  if (length(y) != length(y_pred)) {
    stop("`y` and `y_pred` must have the same length.")
  }

  # common residuals
  err <- y_pred - y

  # --- R2_model ---
  sst <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  sse <- sum((y - y_pred)^2, na.rm = TRUE)
  R2_model <- if (sst > 0) 1 - (sse / sst) else NA_real_

  # --- RMSE / MAE
  RMSE <- sqrt(mean(err^2, na.rm = na_rm))
  MAE  <- mean(abs(err), na.rm = na_rm)

  # --- CV / median
  CV = sd(y_pred, na.rm = na_rm) / mean(y_pred, na.rm = na_rm)
  Median = median(y_pred, na.rm = na_rm)

  # --- handle denominators for relative metrics ---
  denom <- y
  if (rel_zero == "omit") {
    keep <- is.finite(denom) & denom != 0
    rel_err <- (err[keep]) / (denom[keep])
  } else if (rel_zero == "epsilon") {
    rel_err <- err / (denom + epsilon)
  } else { # "none"
    rel_err <- err / denom
  }

  # remove NA/NaN/Inf if requested
  if (na_rm) {
    rel_err <- rel_err[is.finite(rel_err)]
  }

  # if after filtering nothing remains, set relative metrics to NA
  if (length(rel_err) == 0L) {
    RMSRE <- NA_real_
    MARE  <- NA_real_
    MedRelBias <- NA_real_
  } else {
    RMSRE <- sqrt(mean(rel_err^2))
    MARE  <- mean(abs(rel_err))
    MedRelBias <- stats::median(rel_err)
  }

  # --- MedBias (signed) ---
  MedBias <- stats::median(err, na.rm = na_rm)

  # scale to % if requested
  mult <- if (percent) 100 else 1
  RMSRE <- RMSRE * mult
  MARE  <- MARE  * mult
  MedRelBias <- MedRelBias * mult

  # assemble result
  out <- data.frame(
    # R2   = R2_model,
    RMSE       = RMSE,
    MAE        = MAE,
    RMSRE      = RMSRE,
    MARE       = MARE,
    MedBias    = MedBias,
    MedRelBias = MedRelBias,
    Median = Median,
    CV = CV,
    row.names = NULL
  )
  return( out )
}

