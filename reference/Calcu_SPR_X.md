# Calculate F%SPR: the fishing mortality that reduces the spawning biomass per recruit (SBR) to X % of those without fishing (SPR0).

*Calcu_SPR_X* calculates F%SPR (the fishing mortality that reduce the
spawning biomass per recruit (SBR) to X % of those without fishing
(SPR0)).

## Usage

``` r
Calcu_SPR_X(x, M, Sel, w, g, A = 6, F_max = 2)
```

## Arguments

- x:

  The proportion (X) to calculate F% SPR. This parameter should be
  provided as a proportion (i.e. within the range from 0 to 1). Defalut
  is 0.2.

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

  This value should be large enough so that the true F%SPR value be
  included (otherwise a warning message appears). Defalut is 1.

## Value

A value of F% SPR.
