# TMB::MakeADFunに必要な引数から固定効果とランダム効果を推定する関数

TMB::MakeADFunに必要な引数から固定効果とランダム効果を推定する関数

## Usage

``` r
est_mixed(
  tmbdata,
  par_init,
  map,
  random = "U",
  cpp.file.name = "sam2",
  silent = TRUE,
  bias_correct = TRUE
)
```

## Arguments

- tmbdata:

  tmbdata

- par_init:

  パラメータの初期値

- map:

  推定しない（初期値で固定する）パラメータ

- random:

  ランダム効果で推定するパラメータ

- cpp.file.name:

  cppファイルの名前（初期設定："sam"）

- bias_correct:

  平均値のバイアス補正をするかどうか
