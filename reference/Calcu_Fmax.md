# Calculate the fishing mortality that maximizes the yield per recruit.

*Calcu_Fmax* calculates the fishing mortality that maximizes the yield
per recruit.

## Usage

``` r
Calcu_Fmax(M, Sel, w, g, A = 6, F_max = 10, method = "Baranov")
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

- F_max:

  Maximum value of F values to be searched for determining F%SPR. This
  value should be large enough so that the true F%SPR value be included.
  Default is 10

- method:

  Method of the fishing equation. Either "Baranov" or "Pope" is allowed.
  Default is "Baranov".

## Value

A value of the F that maximizes the yield per recruit.
