# 条件付き負の対数尤度を成分別に集計する

SAMの推定結果から、資源尾数過程、漁獲死亡係数過程、および観測系列ごとの
条件付き負の対数尤度を取得します。

## Usage

``` r
get_cond_nll(samres, index_name = NULL)
```

## Arguments

- samres:

  \[sam()\] が返す \`sam\` オブジェクト。

- index_name:

  資源量指標の表示名。\`NULL\` の場合は \`"Index_1"\`, \`"Index_2"\`,
  ... を使用します。漁獲量系列を除く
  観測系列数と同じ長さの文字ベクトルを指定してください。

## Value

次の列を持つ\[tibble\]\[tibble::tibble\]。

- type:

  尤度成分または観測系列を表すfactor。

- nll:

  条件付き負の対数尤度。

## Details

この関数が返す値は、推定されたランダム効果に条件付けた負の対数尤度です。
TMBのLaplace近似によってランダム効果を積分した周辺負の対数尤度では
ありません。また、初期状態の密度、体重・成熟モデル、ランダム効果
\`logB\` の密度、およびペナルティ項は戻り値に含まれません。

C++モデルが返す観測ごとの \`ans_obs\` をfleetごとに合計します。 fleet
1を年齢別漁獲尾数、fleet 2以降を資源量指標として扱います。

## Examples

``` r
if (FALSE) { # \dontrun{
data("samres_example", package = "frasam")
get_cond_nll(samres)
} # }
```
