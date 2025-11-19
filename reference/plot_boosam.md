# ブートストラップについてプロットする関数

ブートストラップについてプロットする関数

## Usage

``` r
plot_boosam(
  samres,
  boores,
  CI = 0.95,
  what.plot = c("biomass", "SSB", "Recruitment", "U"),
  draw_deltaCI = FALSE,
  scenario_name = c("Estiamte", "Bootstrap"),
  alpha = 0.4,
  size = 1,
  base_size = 14,
  log_scale = FALSE,
  legend_name = "Scenario",
  legend_nrow = 1,
  legend_position = "top",
  years = NULL,
  ncol = 2
)
```

## Arguments

- samres:

  [`sam()`](sam.md)の結果オブジェクト

- boores:

  `boo_sam`の結果オブジェクト

- draw_deltaCI:

  デルタ法の信頼区間を描くか（デフォルト:FALSE）

- scenario_name:

  推定値の結果、ブートストラップの結果の凡例に使う名前
