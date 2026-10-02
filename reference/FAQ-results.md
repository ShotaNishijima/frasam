# Example assessment results from the FAQ

Assessment results generated using the examples in `vignettes/FAQ.Rmd`
and the packaged `dat_ex` data.

## Usage

``` r
res_vpa

res_bh
```

## Format

Lists containing assessment inputs, estimated population numbers,
fishing mortality, biomass, spawning biomass, and model diagnostics.

An object of class `vpa` of length 28.

An object of class `sam` of length 34.

## Details

`res_vpa` is the tuned VPA result from the `set-qinit` example, with the
final-year recruitment indices set to missing. `res_bh` is the SAM
result with Beverton-Holt recruitment from the `use-do.call` example.
Saved TMB external pointers cannot be used in a new R session; refit
with `do.call(sam, res_bh$input)` when a live TMB objective is required,
after running
[`use_sam_tmb()`](https://shotanishijima.github.io/frasam/reference/use_sam_tmb.md).

## Examples

``` r
data("res_vpa", package = "frasam")
data("res_bh", package = "frasam")
dim(res_vpa$naa)
#> [1]  7 40
dim(res_bh$naa)
#> [1]  7 40
```
