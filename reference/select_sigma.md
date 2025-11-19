# 観測誤差やプロセス誤差のステップ形式のモデル選択（一つの変数について）

`devide_sigma`関数を順々に実行し、AIC規準で最適なモデルを探索

## Usage

``` r
select_sigma(
  samres,
  var = c("varC", "varF", "varN", "index.key")[1],
  X = NULL,
  stopAIC = TRUE
)
```

## Arguments

- samres:

  sam object

- var:

  分ける誤差の指定.
  "varC"（年齢別漁獲尾数の観測誤差）,"varF"（Fのプロセス誤差）,"varN"（Nのプロセス誤差）,"index.key"（指標の観測誤差）のいずれかを選択．

- X:

  境目を入れる場所の候補

- stopAIC:

  AICが小さくならなかった時点で計算をやめるか（default: TRUE)
