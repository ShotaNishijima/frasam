# Decompose biomass changes from a VPA or SAM result

Extracts numbers at age, fishing mortality, weight, natural mortality,
maturity, and the plus-group setting from a fitted \`vpa\` or \`sam\`
object, then calls \[decompose_biomass_effects()\].

## Usage

``` r
decompose_biomass_factors(
  result,
  target = c("biomass", "ssb"),
  recruitment_age_row = 1L,
  zero_tol = 1e-12
)
```

## Arguments

- result:

  A fitted object of class \`vpa\` or \`sam\`.

- target:

  Quantity to decompose: \`"biomass"\` or \`"ssb"\`.

- recruitment_age_row:

  Row containing recruitment.

- zero_tol:

  Tolerance for zero denominators.

## Value

The result of \[decompose_biomass_effects()\], including age-specific,
age-aggregated, and percentage contributions. Percentages use the total
biomass or SSB in the preceding year as denominator.
