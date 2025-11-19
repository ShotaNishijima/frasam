# 観測誤差やプロセス誤差のステップ形式のモデル選択（複数の変数について）

`devide_sigma`関数を順々に実行し、AIC規準で最適なモデルを探索

## Usage

``` r
select_sigma_grid(
  samres,
  grid = expand.grid(var = c("varC", "varF", "varN", "index.key")[1:2], X = 1:2),
  stopAIC = TRUE,
  check_converge = FALSE,
  SEmax = 10
)
```

## Arguments

- samres:

  sam object

- grid:

  'var'と'X'からなるdata.frame

- stopAIC:

  AICが小さくならなかった時点で計算をやめるか（default: TRUE)
