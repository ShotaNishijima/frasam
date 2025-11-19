# utilities for extracting statistics

utilities for extracting statistics

## Usage

``` r
get_pm(
  res_sam,
  res_future = NULL,
  waa_catch = res_sam$input$dat$waa,
  last_year = 2023,
  waa_biom = res_sam$input$dat$waa,
  year_biol = 2016:2023,
  year_Fcur = 2021:2023,
  perSPR = c(30, 40, 50, 60, 70),
  start_year = 1970
)
```

## Arguments

- res_sam:

  result of sam function

- res_future:

  result of future projection

- waa_catch:

  weight at age for calculating total catch weight. This matrix should
  be given because SAM doesn't use this biological parameters

- last_year:

  last year of stock assessment

- waa_biom:

  weight at age for calculating biomass weight. If this is not given,
  sam_res\$input\$dat\$waa is automatically used

- year_biol:

  year range to derive biological parameters. Biological parameters are
  averaged by age during the given period

- year_Fcur:

  year range to derive F at age considred as F current. Fs at age are
  averaged by age during the given period

- perSPR:

  percent SPR (%) for calculating F%SPR

## Examples

``` r
if (FALSE) { # \dontrun{
res_pm_all <- purrr::map_dfr(basecase_list,
                  function(x) get_pm(x, NULL, waa_catch=x$input$dat$waa), .id="id") %>%
                  pivot_wider(values_from=value, names_from=stat)
} # }
```
