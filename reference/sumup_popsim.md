# Self-test, Cross-testの結果をまとめるための関数

Self-test, Cross-testの結果をまとめるための関数

## Usage

``` r
sumup_popsim(res_true, fit2PS, percent = FALSE, CI = 0.95, ...)
```

## Arguments

- res_true:

  真の推定値をもつSAMまたはVPAのオブジェクト

- fit2PS:

  [`fit2PSdata()`](fit2PSdata.md)
  で得られる、疑似データにフィットさせた結果オブジェクト

- percent:

  Logical; if `TRUE` (default), express relative metrics (RMSRE, MARE,
  MedRelBias) in percentage (i.e., multiplied by 100). If `FALSE`, leave
  as proportions.

## Value

A list with the following components:

- `summary`: A data.frame containing performance metrics. See
  \[calc_metrics()\] for details on the definitions of RMSE, MAE,
  R2_model, etc.

- `all`: すべてのRunについての結果のデータフレーム
