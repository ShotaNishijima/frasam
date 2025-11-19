# Calculate the stock biomass at the maximum sustainable yield (Bmsy).

*Calcu_Bmsy* calculates the stock biomass at the maximum sustainable
yield (Fmsy) for the Bevertoh-Holt stock-recruitment relationship.

## Usage

``` r
Calcu_Bmsy(Fmsy, M, Sel, w, g, A = 6, alpha, beta, method_SR)
```

## Arguments

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

- alpha:

  Estimated Parameter alpha in the Beverton-Holt stock recruitment
  relationship.

- beta:

  Estimated Parameter beta in the Beverton-Holt stock recruitment
  relationship.

## Value

A value of the stock biomass at the maximum sustainable yield (Bmsy).
