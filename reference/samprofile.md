# SAMで推定されたパラメータに対してプロファイル尤度を計算する

SAMで推定されたパラメータに対してプロファイル尤度を計算する

## Usage

``` r
samprofile(
  samres,
  param_name,
  which_param = 1,
  param_range = NULL,
  length = 50
)
```

## Arguments

- samres:

  sam object

- param_name:

  変化させるパラメータの名前

- which_param:

  変化させるパラメータの位置（複数のパラメータが同じ名前を持つときに使用）
