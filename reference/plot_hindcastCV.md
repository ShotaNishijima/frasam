# レトロの結果を使って資源量指標値に対するhindcast cross validationをプロットする関数

レトロの結果を使って資源量指標値に対するhindcast cross
validationをプロットする関数

## Usage

``` r
plot_hindcastCV(
  samres,
  retrores,
  h = 1,
  log = FALSE,
  index_name = NULL,
  show_mase = TRUE,
  mase_position = "upperright",
  use_index = NULL,
  years = NULL
)
```

## Arguments

- samres:

  SAM object

- retrores:

  `retro_sam(res,...)`で実地されたレトロ解析の結果オブジェクト

- log:

  MASEを計算するときに、Indexの観測値と予測値に対してlogを取るかどうか（デフォルトはFALSE）

- index_name:

  各Indexの名前ベクトル、Indexの数だけ必要

- show_mase:

  MASEの結果を載せるかどうか

- use_index:

  特定のIndexを使う場合、`use_index = c(1,3)`のように指定する

- years:

  プロットする期間を指定する場合、`years = 2015:2024`のように指定する
