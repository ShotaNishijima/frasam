# Calculate yield per recruit with a given fishing mortality value.

*Calcu_YPR* calculates yield per recruit from the statistics derived
from the *Calcu_stat* function and input data with a given fishing
mortality value.

## Usage

``` r
Calcu_YPR(Fish_mort, M, Sel, w, g, A = 6, method = "Baranov")
```

## Arguments

- Fish_mort:

  Fishing mortality to calculate yield per recruit. Should be provided
  with a numeric value.

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
  Default is "Baranov".

## Value

A value of yield per recruit with the given fishing mortality.
