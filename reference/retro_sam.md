# レトロスペクティブ解析を実施する

retrospective forecasting is also possible

## Usage

``` r
retro_sam(
  res,
  n = 5,
  stat = "mean",
  b.fix = TRUE,
  remove_short_index = -1,
  map_add = NULL,
  p0_retro_list = NULL
)
```

## Arguments

- res:

  SAM object

- n:

  the number of peels
