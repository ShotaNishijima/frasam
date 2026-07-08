# 条件付き負の対数尤度を成分別に描画する

\[get_cond_nll()\] の出力を横向きの棒グラフで表示します。

## Usage

``` r
plot_cond_nll(x, show_value = TRUE)
```

## Arguments

- x:

  \[get_cond_nll()\] が返すデータフレーム。

- show_value:

  棒の外側に負の対数尤度を小数第1位まで表示するか。 既定値は \`TRUE\`。

## Value

\`ggplot\` オブジェクト。

## Examples

``` r
if (FALSE) { # \dontrun{
data("samres_example", package = "frasam")
plot_cond_nll(get_cond_nll(samres))
} # }
```
