# Calculate F0.1 value.

*Calcu_F0.1* calculates the fishing mortality rate at which the slope of
the yield-per-recruit curve is 10% of the slope of the curve at its
origin.

## Usage

``` r
Calcu_F0.1(M, Sel, w, g, A = 6, method = "Baranov", F_max = 10)
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
  Default is "Baranov".

- F_max:

  Maximum value of F values to be searched for determining F%SPR. This
  value should be large enough so that the true F%SPR value be included.
  Default is 10

## Value

A value of the F0.1.
