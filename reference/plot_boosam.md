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

  [`sam()`](https://shotanishijima.github.io/frasam/reference/sam.md)の結果オブジェクト

- boores:

  `boo_sam`の結果オブジェクト

- CI:

  Confidence interval width. Set 0 to omit confidence intervals.

- what.plot:

  Statistics to plot.

- draw_deltaCI:

  デルタ法の信頼区間を描くか（デフォルト:FALSE）

- scenario_name:

  推定値の結果、ブートストラップの結果の凡例に使う名前

- alpha:

  Alpha value for confidence interval ribbons.

- size:

  Line width.

- base_size:

  Base font size.

- log_scale:

  If `TRUE`, use a log scale on the y-axis.

- legend_name:

  Legend title.

- legend_nrow:

  Number of rows in the legend.

- legend_position:

  Legend position.

- years:

  Years to include in the plot.

- ncol:

  Number of columns in the facet plot.
