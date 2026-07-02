# Check convergence diagnostics for a SAM fit

`check_fit_sam` summarizes basic convergence diagnostics for a
[`sam()`](https://shotanishijima.github.io/frasam/reference/sam.md)
result object. The function does not change the fitted object.

## Usage

``` r
check_fit_sam(
  res,
  gradient_tol = 0.01,
  se_max = Inf,
  sigma_range = c(1e-04, 10),
  rho_range = c(1e-04, 1 - 1e-04),
  par_abs_max = Inf,
  boundary_tol = 1e-04,
  verbose = TRUE
)
```

## Arguments

- res:

  A `sam` result object.

- gradient_tol:

  Maximum allowed absolute fixed-effect gradient.

- se_max:

  Maximum allowed fixed-effect standard error. Set `Inf` to skip this
  check.

- sigma_range:

  Allowed range for reported standard deviations.

- rho_range:

  Allowed range for reported rho values.

- par_abs_max:

  Maximum allowed absolute fixed-effect estimate on the internal scale.
  Set `Inf` to skip this check.

- boundary_tol:

  Tolerance for detecting parameters close to finite lower or upper
  bounds.

- verbose:

  If `TRUE`, print the diagnostic table.

## Value

A list with `ok`, `checks`, `fixed`, `sigma`, and `boundary` elements.
