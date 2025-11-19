# 観測誤差とプロセス誤差をどこかの年齢間で分けて推定し直す関数

観測誤差とプロセス誤差をどこかの年齢間で分けて推定し直す関数

## Usage

``` r
divide_sigma(
  samres,
  var = c("varC", "varF", "varN", "index.key")[1],
  which = (1:6)[1]
)
```

## Arguments

- samres:

  sam object

- var:

  分ける誤差の指定.
  "varC"（年齢別漁獲尾数の観測誤差）,"varF"（Fのプロセス誤差）,"varN"（Nのプロセス誤差）,"index.key"（指標の観測誤差）のいずれかを選択．

- which:

  どこで区切りをいれるか.
  1だと0歳と1歳以上で分け（加入が0歳の場合）、指標値の場合は1本目と2本目以降の観測誤差を分ける
