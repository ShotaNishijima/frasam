# Calculate the fishing mortality that maximizes the maximum sustainable yield (Fmsy) with fixed steepness parameters (h) in the Bevertoh-Holt stock-recruitment relationship.

*Calcu_Fmsy* calculates the fishing mortality that maximizes the maximum
sustainable yield (Fmsy) with fixed steepness parameters (h).

## Usage

``` r
Calcu_Fmsy(M, Sel, w, g, A = 6, method, F_max = 10, alpha, beta, method_SR)
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

- method:

  Method of the fishing equation. Either "Baranov" or "Pope" is allowed.

- F_max:

  Maximum value of F values to be searched for determining Fmsy.

- alpha:

  Estimated Parameter alpha in the Beverton-Holt stock recruitment
  relationship.

- beta:

  Estimated Parameter beta in the Beverton-Holt stock recruitment
  relationship.

## Value

A value of the F that maximizes the yield per recruit.
