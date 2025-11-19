# レトロの結果から各Indexに対する予測値を抽出して、Mean Absolute Scaled Errorを計算する関数

レトロの結果から各Indexに対する予測値を抽出して、Mean Absolute Scaled
Errorを計算する関数

## Usage

``` r
calc_mase(samres, retrores, h = 1, log = FALSE, index_name = NULL)
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
