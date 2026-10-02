# SAMによる資源計算を実施する

SAMによる資源計算を実施する

## Usage

``` r
sam(
  dat,
  last.catch.zero = FALSE,
  cpp.file.name = "sam2",
  tmb.run = FALSE,
  abund = c("B"),
  rec.age = 0,
  min.age = rep(0, length(abund)),
  max.age = rep(0, length(abund)),
  plus.group = TRUE,
  alpha = 1,
  b.est = FALSE,
  b.fix = rep(NA, length(abund)),
  varC = 0,
  varN = 0,
  varF = 0,
  varN.fix = NULL,
  est.method = NULL,
  SR = "BH",
  AR = 0,
  rho.mode = 2,
  q.init = NULL,
  sdFsta.init = NULL,
  sdLogN.init = NULL,
  sdLogObs.init = NULL,
  rho.init = NULL,
  a.init = NULL,
  b.init = NULL,
  ref.year = 1:5,
  bias.correct = TRUE,
  bias.correct.sd = FALSE,
  get.random.vcov = FALSE,
  silent = FALSE,
  remove.Fprocess.year = NULL,
  RW.Forder = 0,
  map.add = NULL,
  p0.list = NULL,
  scale = 1000,
  gamma = 10,
  sel.def = "max",
  use.index = NULL,
  upper = NULL,
  lower = NULL,
  index.key = NULL,
  index.b.key = NULL,
  b_random = FALSE,
  b_range = NULL,
  lambda = 0,
  FreeADFun = FALSE,
  add_random = NULL,
  lambda_Mesnil = 0,
  tmbdata = NULL,
  map = NULL,
  model_wm = c(FALSE, FALSE),
  w0_factor = c("none", "SSB")[1],
  weight_factor = c("none", "N_total")[1],
  family_w = c("lognormal", "gamma")[2],
  maturity_factor = c("none", "cohort_plus")[1],
  scale_number = 1000,
  weight_weight = NULL,
  maturity_weight = NULL,
  g_fix = NULL,
  CV_w_fix = NULL,
  w_link = "log",
  sep_omicron = TRUE,
  growth_regime = NULL,
  catch_prop = NULL,
  no_est = FALSE,
  getJointPrecision = FALSE,
  loopnum = 2,
  obj_overwrite = NULL,
  ignore.parm.uncertainty = FALSE,
  change_plusgroup = FALSE
)
```

## Arguments

- dat:

  samに使用するdataでrvpaと同じフォーマットで利用可能

- last.catch.zero:

  最終年の漁獲量がない場合。デフォルトはFALSE（最終年の漁獲量が利用できて用いる）

- cpp.file.name:

  推定に用いるcppファイル。デフォルトは最新版の"sam2"

- abund:

  Indexの種類。用いるIndexの長さのベクトル。"B":
  資源量、"SSB"：親魚資源量、"N"：尾数、"Bs"：資源量×全体のFから計算した選択率、"Bf"：資源量×指標ごとの選択率（Fを`catch_prop`で調整）。合計する年齢の幅については`min.age`と`max.age`。オメガ法との対応はDetailsを参照。

- rec.age:

  加入年齢 (default:
  0)。1以上の場合でも、min.ageやmax.ageは0歳始まりとしての年齢を入れること。

- min.age:

  Indexの最低年齢
  ([`frasyr::vpa()`](https://rdrr.io/pkg/frasyr/man/vpa.html)と同じで最小の年齢を0とする)
  用いるIndexの長さのベクトル

- max.age:

  Indexの最高年齢
  ([`frasyr::vpa()`](https://rdrr.io/pkg/frasyr/man/vpa.html)と同じで最小の年齢を0とする)
  用いるIndexの長さのベクトル

- plus.group:

  プラスグループを考慮する（TRUE, デフォルト）、考慮しない（FALSE）

- alpha:

  最高年齢-1歳へのFに対する最高年齢のFの比

- b.est:

  Indexと資源量の間の非線形関係を考慮しない（FALSE,
  デフォルト）、考慮する(TRUE)

- b.fix:

  b.est=TRUEの場合、非線形パラメータbを推定するか（NA）、固定するか（固定する値）
  ??

- varC:

  CAAの観測誤差

- varN:

  ??

- varF:

  ??

- est.method:

  Deprecated. Use `index.key` instead. If `"ls"` is supplied and
  `index.key = NULL`, all index observation-error sigmas are shared.

- SR:

  再生産関係："RW"(Random walk), "BH", "RI", "HS", "Mesnil", or "Const"

- AR:

  再生産関係の残差の自己相関パラメータを推定する（１）、推定しない（0,
  デフォルト） ?? これでいい？

- rho.mode:

  Fの多変量Random
  walkの非対角成分の相関係数のタイプ。１ならすべての年齢間で相関係数1,
  0なら完全にランダム,
  2は任意の年齢間で共通のrhoを推定、3は年齢i,jの相関を\\rho^\|i-j\|\\で推定

- q.init:

  パラメータqの初期値。NULLの場合にはexp(-5)が用いられる。

- sdFsta.init:

  Fのランダムウォークのσの初期値。NULLの場合（デフォルト）にはlog_sigma=-0.693147が用いられる。初期値を与える場合には、ノーマルスケールでの値を与える。

- sdLogN.init:

  Nのプロセス誤差のσの初期値。NULLの場合（デフォルト）にはlog_sigma=0.35が用いられる（←コード上はこうなっているがこれで良い？）。初期値を与える場合には、ノーマルスケールでの値を与える。

- sdLogObs.init:

  Indexの観察誤差の初期値。NULLの場合（デフォルト）にはlog_sigma=-0.356675が用いられる。初期値を与える場合には、ノーマルスケールでの値を与える。

- rho.init:

  rhoの初期値。NULLの場合（デフォルト）には0が用いられる。

- a.init:

  再生産関係パラメータaの初期値。NULLの場合（デフォルト）、log_a=8
  (SR="Const")またはlog_a=-4（SR="Const"以外）

- b.init:

  再生産関係パラメータbの初期値。NULLの場合（デフォルト）、log_b=7
  (SR="HS", "Mesnil", "BHS") またはlog_b=-8（それ以外）

- ref.year:

  管理基準値を計算するときの参照年。最終年から何年分さかのぼるか。デフォルトは1:5（最新年からさかのぼって５年分）。

- bias.correct:

  固定効果パラメータ？？のバイアスを補正する（TRUE：デフォルト）、補正しない（FALSE）

- bias.correct.sd:

  ランダム効果のSDパラメータ？？のバイアスを補正する（TRUE）、補正しない（FALSE：デフォルト）

- get.random.vcov:

  ランダム効果の分散共分散行列を推定する（TRUE、時間かかります）、推定しない（FALSE：デフォルト）

- silent:

  MakeADfunのときの標準出力あり（TRUE: デフォルト）、なし（FALSE）

- remove.Fprocess.year:

  Years to remove from the F process model.

- RW.Forder:

  Order of the random-walk process for fishing mortality.

- map.add:

  Additional TMB map settings.

- p0.list:

  Initial parameter list.

- scale:

  資源量のスケーリングファクター。資源量はscaleで割った値となる

- gamma:

  Gamma parameter setting.

- sel.def:

  選択率の定義。"max"（デフォルト）は最大値で、"mean"は全年齢の平均値で、"maxage"は最高年齢の値でFを割る。"Bf"では`F * catch_prop`を同様に正規化する。

- use.index:

  この仕様は設定ミスを引き起こしやすいので廃止しました。`sam()`を使用するまえに、Indexのデータを使用するもののみにsubsetするようにしてください

- upper:

  推定パラメータの上限値。固定効果の数のLengthを持つ必要あり。NULL（デフォルト）の場合Inf

- lower:

  推定パラメータの下限値。固定効果の数のLengthを持つ必要あり。NULL（デフォルト）の場合-Inf

- index.key:

  Index間の観測誤差sigmaの制約。NULLの場合はIndexごとに別々、rep(0,
  length(abund))の場合は全Indexで共通。

- index.b.key:

  Indexのbの制約 ?? どうやって使う？？

- b_random:

  再生産パラメータbをランダム効果として推定するかどうか??

- b_range:

  再生産パラメータbの値の範囲

- lambda:

  ??

- FreeADFun:

  [`TMB::FreeADFun`](https://rdrr.io/pkg/TMB/man/FreeADFun.html)を使う場合、TRUEにする。See
  [`?TMB::FreeADFun`](https://rdrr.io/pkg/TMB/man/FreeADFun.html).

- add_random:

  Additional random effects.

- lambda_Mesnil:

  SRを"Mesnil"にした場合のラムダの値??

- tmbdata:

  TMB data list.

- map:

  TMB parameter map.

- model_wm:

  weightとmaturityの成長をモデリングするかどうか

- w0_factor:

  Factor setting for initial weight.

- weight_factor:

  Factor setting for weight-at-age.

- family_w:

  Error distribution for weight observations.

- maturity_factor:

  Factor setting for maturity-at-age.

- scale_number:

  尾数のスケーリングファクター。尾数はscaleで割った値となる

- weight_weight:

  Weight assigned to weight-at-age likelihood components.

- maturity_weight:

  Weight assigned to maturity-at-age likelihood components.

- g_fix:

  Fixed value setting for growth parameter g.

- CV_w_fix:

  Fixed value setting for weight coefficient of variation.

- w_link:

  Link function for weight model.

- sep_omicron:

  Whether to separate omicron parameters.

- growth_regime:

  Growth regime setting.

- catch_prop:

  `abund = "Bf"`のときにFに掛ける年齢・年・指標別の重み。
  通常は対象fleetのcatch at
  ageの全体に対する比率。オメガ法ではomegaを指定できる。
  次元は`c(nrow(dat$waa), ncol(dat$waa), nrow(dat$index) + 1)`で、
  第3次元の1番目は漁獲データ、`i + 1`番目は`i`番目のIndexに対応する。
  `NULL`の場合はすべて1。

- no_est:

  推定しない(TRUE)、パラメータ推定する（FALSE, デフォルト）

- getJointPrecision:

  JointPrecision matrixを計算しない（FALSE,
  デフォルト）、計算する（TRUE)

- loopnum:

  最適化を繰り返す回数．デフォルトは2

- change_plusgroup:

  プラスグループが途中で変わるか

- varNfix:

  ??

## Details

`abund = "Bf"`は、`catch_prop`の対象Indexのスライスにomegaを指定し、
`sel.def = "mean"`とすることで、`frasyr::vpa(abund = "Bo")`
と実質的に同じ指標への当てはめ方に対応する。
VPAの"Bo"は`F * omega`を年ごとに合計が1となるよう正規化する。
SAMの"Bf"で平均値による正規化を使うと、その重みはVPAの年齢区分数倍となる。
年齢区分数が一定でqを推定する場合、この定数倍はqに吸収される
（資源量のスケーリングもqに影響する）。VPAとSAMの推定結果自体が一致するという意味ではない。
同じN、F、体重、omega、対象年齢を用いた場合、SAMとVPAのscaleをそれぞれ
`scale_sam`、`scale_vpa`、全年齢区分数を`K`とすると、
`A_Bf = (K * scale_vpa / scale_sam) * A_Bo`となる。
指標の関係が`I = q * A^b`の場合、Boの尺度への換算は
`q_on_Bo_scale = q_Bf * (K * scale_vpa / scale_sam)^b`。
同じscaleで`b = 1`なら`q_on_Bo_scale = K * q_Bf`となる。
omega全体の定数倍は選択率の正規化で打ち消されるため、omegaを割るだけではqの尺度は揃わない。
この換算は正規化・単位の違いを補正するもので、別々に推定したVPAとSAMのqの一致を保証しない。
`sel.def = "max"`では正規化の比率が年ごとに変わり得るため、
一定のqで吸収できるとは限らない。`sel.def`は他の"Bs"・"Bf"にも共通して適用される。
詳細は[FAQのオメガ法の項目](https://shotanishijima.github.io/frasam/articles/FAQ.html#omega-method)を参照。
