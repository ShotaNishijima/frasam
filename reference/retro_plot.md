# レトロスペクティブ解析の結果をプロットする

レトロスペクティブ解析の結果をプロットする

## Usage

``` r
retro_plot(
  res,
  retro_res,
  start_year = NULL,
  scale = 1000,
  forecast = FALSE,
  base_size = 14,
  plot_mohn = TRUE,
  mohn_position = "upperleft"
)
```

## Arguments

- res:

  sam object

- retro_res:

  `retro_sam`の結果

- start_year:

  プロットを開始する年
