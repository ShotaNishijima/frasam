# SAMの固定効果だけを推定する関数

SAMの固定効果だけを推定する関数

## Usage

``` r
est_fixed(
  tmbdata,
  par_init,
  map,
  random = "U",
  cpp.file.name = "sam2",
  silent = TRUE,
  ...
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
