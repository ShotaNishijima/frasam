# Calculate the relative number at age at the equilibrium.

*Calcu_N* calculates the relative number at age at the equilibrium (N)
with a given fishing mortality value.

## Usage

``` r
Calcu_N(Fish_mort, M, Sel, A = 6)
```

## Arguments

- Fish_mort:

  Fishing mortality to calculate SBR. Should be provided with a numeric
  value.

- M:

  A vector of natural mortality rate at age a. The length should be A+1.

- Sel:

  A vector of selectivity at age a. The length should be A+1.

- A:

  Plus group age. Default is 6.
