# Decompose biomass or spawning-stock-biomass changes by age

Decompose biomass or spawning-stock-biomass changes by age

## Usage

``` r
decompose_biomass_effects(
  naa,
  waa,
  faa,
  M,
  maa = NULL,
  plus_group = TRUE,
  recruitment_age_row = 1L,
  zero_tol = 1e-12
)
```

## Arguments

- naa:

  Matrix of numbers at age. Rows are ages, columns are years.

- waa:

  Matrix of body weights at age. Rows are ages, columns are years.

- faa:

  Matrix of fishing mortality at age. Rows are ages, columns are years.

- M:

  Matrix or vector of natural mortality at age.

- maa:

  Optional matrix of maturity at age. If NULL, biomass is decomposed. If
  supplied, spawning-stock biomass is decomposed.

- plus_group:

  Logical. If TRUE, the last row is treated as a plus group. If higher
  age rows are missing for some years, ages above the common plus-group
  row are collapsed before decomposition.

- recruitment_age_row:

  Row index for recruitment. Default is 1, assumed to be age 0.

- zero_tol:

  Tolerance for zero denominators.

## Value

A list of matrices with the same age-year dimensions used for
decomposition. When the plus group is harmonized, the matrices have the
collapsed common plus-group age rows.
