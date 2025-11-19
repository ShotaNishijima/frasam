# Estimate the Beverton-Holt and hockey-stick stock recruitment relationships from data of spawning stock biomass (SBy) and recruits (Ry) with a give steepness (h) value.

*Est_SR* estimates the coefficients (i.e. alpha and beta) of the BH
relationship from data of spawning stock biomass (SBy) and recruits
(Ry).

## Usage

``` r
Est_SR(
  SBy,
  Ry,
  h,
  h_fix = T,
  inits = c(1, 1),
  ref.year = c(1970:2019),
  SB0_max = max(SBy) * 100,
  M,
  w,
  g,
  A = 6,
  method
)
```

## Arguments

- SBy:

  A vector of spawning biomass at year. The vector should be named with
  years. The length should be the same as that of *ref.year*

- Ry:

  A vector of recruits at year. The vector should be named with years.
  The length should be the same as that of *ref.year*

- h:

  A value of steepness to be fixed.

- h_fix:

  A logical value deciding whether the steepness is fixed. Default is
  TRUE. FALSE is allowed only for the BH relationship at this stage.

- inits:

  A vector of initial values of the parameters search (alpha, beta) for
  BH. Although default is c(1, 1), searching for the suitable sets is
  recommended.

- ref.year:

  A vector of years to be referred for the estimation. Default is
  c(1970:2019)

- SB0_max:

  Maximum value of SB0 values to be searched for parameter estimation.
  Default is max(SBy)\*100.

- M:

  A vector of natural mortality rate at age a. The length should be A+1.

- w:

  A vector of weight at age a. The length should be A+1.

- g:

  A vector of maturity at age a. The length should be A+1.

- A:

  Plus group age. Default is 6.

- method:

  The type of the stock-recruitment relationship. At this stage, "BH"
  (Beverton-Holt) or "HS" (Hockey-stick) is allowed.

- Fish_mort:

  Fishing mortality to calculate SBR. Should be provided with a numeric
  value.

## Value

Values of alpha and beta coefficients for the BH relationships or a
value of SB0 for the HS relationships. As for the BH, when h is not
fixed (i.e. h_fix = T), the output from the optim() function is shown.
