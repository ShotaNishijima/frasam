# Popsimの結果ををプロットする関数

Popsimの結果ををプロットする関数

## Usage

``` r
plot_popsim(
  res_true,
  fit2PS,
  CI = 0.95,
  what.plot = c("biomass", "SSB", "Recruitment", "U"),
  Age = NULL,
  scenario_name = c("True", "Simulation"),
  alpha = 0.4,
  size = 1,
  base_size = 14,
  log_scale = FALSE,
  legend_name = "Scenario",
  legend_nrow = 1,
  legend_position = "top",
  years = NULL,
  ncol = 2,
  sim_colour = "blue"
)
```

## Arguments

- res_true:

  真の推定値をもつSAMまたはVPAのオブジェクト

- fit2PS:

  [`fit2PSdata()`](https://shotanishijima.github.io/frasam/reference/fit2PSdata.md)
  で得られる、疑似データにフィットさせた結果オブジェクト
