# SAM or VPAの結果を描くグラフ

SAM or VPAの結果を描くグラフ

## Usage

``` r
plot_samvpa(
  vpa_sam_list,
  CI = 0.95,
  scenario_name = NULL,
  alpha = 0.4,
  size = 1,
  base_size = 14,
  log_scale = FALSE,
  legend_name = "Scenario",
  legend_nrow = 1,
  legend_position = "top",
  what.plot = c("biomass", "SSB", "Recruitment", "U"),
  years = NULL,
  ncol = 2
)
```

## Arguments

- vpa_sam_list:

  VPAまたはSAMの結果オブジェクトのリスト

- what_plot:

  どの統計量をプロットするか
