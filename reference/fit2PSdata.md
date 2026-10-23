# 生成された疑似データにVPA/SAMを推定し、Self-test/Cross-testを実行

生成された疑似データにVPA/SAMを推定し、Self-test/Cross-testを実行

## Usage

``` r
fit2PSdata(
  res,
  PSdata,
  reest_lambda = FALSE,
  target_retro = NULL,
  n_retro = 5,
  b_fix = TRUE,
  ...
)
```

## Arguments

- res:

  フィットさせるSAM/VPAのオブジェクト

- PSdata:

  `popsim_vpa()`で生成された疑似データ

- reest_lambda:

  ridgeVPAで毎回lambdaを再推定するかどうか．TRUEの場合はlambda探索を実行する

- target_retro:

  ペナルティλの選択に用いるmohn's rhoのパラメータ

- n_retro:

  レトロスペクティブ解析で遡る年数。デフォルトは\`5\`

- b_fix:

  レトロスペクティブ解析内でbを固定するか。デフォルトは\`TRUE\`

- ...:

  frasyr::autocalc_ridgevpaで設定可能な引数
