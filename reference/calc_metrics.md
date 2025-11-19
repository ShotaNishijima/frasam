# Calculate regression/forecast accuracy and bias metrics

Calculates a set of accuracy (error) and bias metrics between observed
values and predictions:

- **R2\\model**: \\1 - (y - y)^2 (y - y)^2\\

- **RMSE**: \\mean((y-y)^2)\\

- **MAE**: \\mean(\|y-y\|)\\

- **RMSRE**: \\mean(((y-y)/y)^2)\\ (optionally in %)

- **MARE**: \\mean(\|y-y\|/y)\\ (optionally in %)

- **MedBias**: \\median(y-y)\\

- **MedRelBias**: \\median((y-y)/y)\\ (optionally in %)

Relative metrics (RMSRE/MARE/MedRelBias) use the observed values `y` in
the denominator. If zeros occur in `y`, control the behavior via
`rel_zero`:

- `"omit"`: drop pairs where \\y=0\\ (default for stability)

- `"epsilon"`: replace \\y\\ by \\y + \\ with small `epsilon`

- `"none"`: no guard; may yield `Inf`/`NaN`

## Usage

``` r
calc_metrics(
  y,
  y_pred,
  na_rm = TRUE,
  percent = FALSE,
  rel_zero = c("omit", "epsilon", "none"),
  epsilon = 1e-08
)
```

## Arguments

- y:

  Numeric vector of observed values.

- y_pred:

  Numeric vector of predicted values (same length as `y`).

- na_rm:

  Logical; if `TRUE` (default), remove `NA`/`NaN`/`Inf` after forming
  needed quantities for each metric.

- percent:

  Logical; if `TRUE` (default), express relative metrics (RMSRE, MARE,
  MedRelBias) in percentage (i.e., multiplied by 100). If `FALSE`, leave
  as proportions.

- rel_zero:

  Character; handling of zeros in `y` for relative metrics. One of
  `"omit"`, `"epsilon"`, `"none"`. See Details.

- epsilon:

  Numeric; small constant added to denominators when
  `rel_zero = "epsilon"`. Default `1e-8`.

## Value

A one-row `data.frame` with columns:
`R2_model, RMSE, MAE, RMSRE, MARE, MedBias, MedRelBias`.

## Details

Compute common regression metrics (R2, RMSE, MAE, RMSRE, MARE, MedBias,
MedRelBias)

- **R2\\model** uses the usual in-sample definition with the sample mean
  of `y`. If `var(y) == 0`, `R2_model` is returned as `NA`.

- Bias-type metrics retain the sign (systematic over-/under-prediction),
  while error-type metrics use absolute value or squares (magnitude
  only).

## Examples

``` r
y      <- c(10, 12, 9, 15, 0, 8)
y_pred <- c(11, 11, 8, 14, 0.2, 9)

# Default: percent = TRUE, rel_zero = "omit"
calc_metrics(y, y_pred)
#>        RMSE       MAE     RMSRE       MARE MedBias  MedRelBias Median        CV
#> 1 0.9165151 0.8666667 0.0993575 0.09722222    -0.4 -0.06666667     10 0.5321906

# Use proportions (not %) and guard zeros with epsilon
calc_metrics(y, y_pred, percent = FALSE, rel_zero = "epsilon", epsilon = 1e-6)
#>        RMSE       MAE    RMSRE     MARE MedBias MedRelBias Median        CV
#> 1 0.9165151 0.8666667 81649.66 33333.41    -0.4 0.01666666     10 0.5321906

# No guarding (may produce Inf/NaN if y contains zeros)
calc_metrics(y, y_pred, rel_zero = "none")
#>        RMSE       MAE     RMSRE       MARE MedBias  MedRelBias Median        CV
#> 1 0.9165151 0.8666667 0.0993575 0.09722222    -0.4 -0.06666667     10 0.5321906
```
