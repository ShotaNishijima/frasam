# Calculate spawning biomass per recruit (SBR) with a given fishing mortality value.

*Calcu_SBR* calculates SBR from the statistics derived from the
*Calcu_stat* function and input data with a given fishing mortality
value.

## Usage

``` r
Calcu_SBR(Fish_mort, M, Sel, w, g, A = 6)
```

## Arguments

- Fish_mort:

  Fishing mortality to calculate SBR. Should be provided with a numeric
  value.

- M:

  A vector of natural mortality rate at age a. The length should be A+1.

- Sel:

  A vector of selectivity at age a. The length should be A+1.

- w:

  A vector of weight at age a. The length should be A+1.

- g:

  A vector of maturity at age a. The length should be A+1.

- A:

  Plus group age. Default is 6.

## Value

A value of spawing biomass par recruit with the given fishing mortality.
