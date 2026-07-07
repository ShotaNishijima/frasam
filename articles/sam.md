# SAMを使った資源量推定

## SAMのインストール

- インストールに必要なパッケージ(devtools)はあらかじめインストールしておいてください

``` r


## インストール (最新ブランチはcmsa11→vignetteとかが追加されたcreate_vignette)
#devtools::install_github("ShotaNishijima/frasam@create_vignette") # frasam
## frasyrも使うのでfrasyrもインストールしてください
#devtools::install_github("ichimomo/frasyr@dev")

# ライブラリの呼び出し
library(frasyr)
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(export_all = FALSE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all()
} else {
  library(frasam)
}
library(tidyverse)
library(patchwork)
# install.packages("lemon")
library(lemon)
```

## SAMを適用するデータセットの作成

VPAを適用するときのデータセットと同じ形式のデータセットが利用できます。
ただし、VPAでは資源量指数が利用できなくても適用できますが、SAMでは資源量指数の利用が必須になります。

### データの読み込み

``` r


caa   <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex1_caa.csv",  row.names=1)
waa   <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex1_waa.csv",  row.names=1)
maa   <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex1_maa.csv",  row.names=1)
index <- read.csv("https://raw.githubusercontent.com/ichimomo/frasyr/dev/data-raw/ex1_index.csv",  row.names=1)
# create pseudo data for caa
set.seed(1)
caa   <- caa * exp(rnorm(length(caa), mean=-0.5 * (0.1)^2, sd=0.1))

# Indexがひとつだと例としてあまり適切でないので、適当な加入Indexを作成する
rec_idx <- seq(1,0.25,length=10)*exp(rnorm(10,0,0.2))
index[2,] <- rec_idx

dat <- data.handler(caa=caa, waa=waa, maa=maa, M=0.5, index=index)
```

## SAMの解析

### 基本的な設定で解析と結果の出力

``` r

# まずはsamのhelpファイルを閲覧してどんなオプションがあるか把握しましょう
#help(sam)

# SAMを使う準備 (samのオプションtmb.run=TRUEでも良いかもだけど、今はうまく動かないみたい？)
use_sam_tmb(TmbFile="sam2",overwrite=FALSE) #compile and activate
#> [1] TRUE

# samの実施
# (*)の部分が設定ポイント
res_sam <- sam(dat,
               last.catch.zero = FALSE,
               cpp.file.name = "sam2",
               tmb.run = FALSE, # TRUEにするとtmbをコンパイルしてロードする。最初の１回だけTRUEにするのが良い？sam2だと動かない
               rec.age = 0,
               plus.group = TRUE,
               alpha = 1, # 最高年齢と最高年齢-1歳のFが同じと仮定（VPAと同じ仮定)
               # est.method = "ml",
               # 資源量指数に関わる設定。資源量指数の数だけベクトルで与える。
               abund = c("SSB", "N"),
               min.age = c(0, 0),
               max.age = c(6, 0),
               b.est = FALSE, # b推定の有無(*)
               b.fix = NA,
               # 分散に関わる設定。どの分散を同じとするか。年齢の数だけ与える (*)
               varC =  c(0,0,0,0,0,0,0),　#とりあえずすべての年齢で共通
               varN = c(0,1,1,1,1,1,1),　#0歳と1歳魚以上で分ける
               varF = c(0,0,0,0,0,0,0),   #とりあえずすべての年齢で共通
               # 1歳魚以上のlog(N)のプロセス誤差を小さい値で固定(分散の値なので、この場合SD=0.01に相当)
               varN.fix = c(NA, 1e-4),
               # 再生産関係推定に関わる設定 (*)
               SR = "BH",
               AR = 0,
               # Fの多変量Random walkの非対角成分の相関係数のタイプ。１ならすべての年齢間で相関係数1, 0なら完全にランダム, 2は任意の年齢間で共通のrhoを推定、3は年齢i,jの相関を\eqn{\rho^|i-j|}で推定 (*)
               rho.mode = 0,
               # 初期値の設定
               q.init = NULL,
               sdFsta.init = NULL,
               sdLogN.init = c(0.2, 0.2),
               sdLogObs.init = NULL,
               rho.init = NULL,
               a.init = NULL,
               b.init = NULL,
               # その他
               # ref.year = 1:5, # 重要そうだけど使われていない
               bias.correct = TRUE,
               bias.correct.sd = FALSE,
               get.random.vcov = FALSE,
               silent = TRUE,
               remove.Fprocess.year = NULL,
               RW.Forder = 0,
               map = NULL, # map (パラメータの推定の有無を個別に調整する) を入れる
               map.add = NULL, # mapを追加する
               p0.list = NULL,
               scale = 1000,
               scale_number = 1000,
               gamma = 10,
               sel.def = "max",
               use.index = NULL,
               upper = NULL,
               lower = NULL,
               index.key = NULL,
               index.b.key = NULL,
               lambda = 0,
               add_random = NULL,
               tmbdata = NULL,
               sep_omicron = TRUE,
               catch_prop = NULL,
               no_est = FALSE,
               getJointPrecision = FALSE,
               loopnum = 2
)

## 推定結果の確認: opt=推定パラメータや目的関数、収束の有無
## - $convergenceが0なら収束
res_sam$opt
#> $par
#>         logQ         logQ logSdLogFsta    logSdLogN  logSdLogObs  logSdLogObs 
#>     9.560474    -7.232696    -1.877590    -1.761564    -3.101510    -1.367962 
#>  logSdLogObs     rec_loga     rec_logb 
#>    -1.082575     2.365818   -20.967589 
#> 
#> $objective
#> [1] -4.304594
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 87
#> 
#> $evaluations
#> function gradient 
#>      116       88 
#> 
#> $message
#> [1] "relative convergence (4)"

## 推定結果の確認: rep=推定パラメータとSE
## - 推定値(Estimate)に対してStd.Errorが異常に大きくないか
## - Maximum gradient componentが十分小さい値か（1e-3以下なら十分OK、多数のパラメータを推定する難しいモデルなら上限0.1以下くらいでも？)
res_sam$rep
#> sdreport(.) result
#>                Estimate   Std. Error
#> logQ           9.560474 8.811785e-02
#> logQ          -7.232696 1.118901e-01
#> logSdLogFsta  -1.877590 1.476098e-01
#> logSdLogN     -1.761564 2.699047e-01
#> logSdLogObs   -3.101510 7.511999e-01
#> logSdLogObs   -1.367962 2.342732e-01
#> logSdLogObs   -1.082575 2.288350e-01
#> rec_loga       2.365818 6.530281e-02
#> rec_logb     -20.967589 4.742059e+04
#> Maximum gradient component: 1.00355e-05

## 推定結果の確認: AIC（単独ではあまり意味がない。モデルを比較するときに）
res_sam$aic
#> [1] 9.390811

## 推定結果の確認: 分散パラメータのsigma。ノーマルスケールで、大きさを確認。
## - すごく小さい値で推定されている場合には推定の意味がないので推定しないなどする
## - どのsigmaを共通にするか？を決める
res_sam[c("sigma","sigma.logC","sigma.logFsta","sigma.logN")]
#> $sigma
#> [1] 0.2546253 0.3387221
#> 
#> $sigma.logC
#> [1] 0.04498123 0.04498123 0.04498123 0.04498123 0.04498123 0.04498123 0.04498123
#> 
#> $sigma.logFsta
#> [1] 0.1529583 0.1529583 0.1529583 0.1529583 0.1529583 0.1529583 0.1529583
#> 
#> $sigma.logN
#> [1] 0.171776 0.010000 0.010000 0.010000 0.010000 0.010000 0.010000

## 固定効果の確認・出力
fixef = out_par(res_sam,filename="FEpar")
knitr::kable(fixef)
```

| FE           |        MLE |           SE | Gradient | Unlinked value |
|:-------------|-----------:|-------------:|---------:|---------------:|
| logQ         |   9.560474 | 8.811780e-02 | -5.4e-06 |   1.419257e+04 |
| logQ         |  -7.232697 | 1.118901e-01 |  1.0e-05 |   7.226000e-04 |
| logSdLogFsta |  -1.877590 | 1.476098e-01 |  1.8e-06 |   1.529583e-01 |
| logSdLogN    |  -1.761564 | 2.699047e-01 |  2.3e-06 |   1.717760e-01 |
| logSdLogObs  |  -3.101510 | 7.511999e-01 | -4.0e-07 |   4.498120e-02 |
| logSdLogObs  |  -1.367962 | 2.342732e-01 | -7.0e-07 |   2.546253e-01 |
| logSdLogObs  |  -1.082575 | 2.288350e-01 | -4.5e-06 |   3.387221e-01 |
| rec_loga     |   2.365818 | 6.530280e-02 |  1.3e-06 |   1.065275e+01 |
| rec_logb     | -20.967589 | 4.742059e+04 |  0.0e+00 |   0.000000e+00 |

``` r


## 結果のプロット
# デルタ法にかかる信頼区間が描かれる
plot_samvpa(res_sam, CI=0.95)
```

![](sam_files/figure-html/initial%20analysis-1.png)

``` r


## 結果の出力
out_sam(res_sam, filename="sam")
```

### VPAや、異なる設定のSAMとの比較

``` r

## VPAもやってみる
res_vpa <- vpa(dat,fc.year=1998:2000,tf.year = 1998:1999,
               term.F="max",stat.tf="mean",Pope=TRUE,tune=TRUE,p.init=0.5, abund=c("SSB","N"), min.age=c(0,0), max.age=c(6,0), sel.update=TRUE)

## VPAとSAMのモデルの比較
gg_graph1 <- plot_vpa(list(VPA=res_vpa, SAM=res_sam))
gg_graph2 <- plot_vpa(list(VPA=res_vpa, SAM=res_sam), what.plot=c("fishing_mortality"))

print(gg_graph1)
```

![](sam_files/figure-html/comparison-1.png)

``` r

print(gg_graph2)
```

![](sam_files/figure-html/comparison-2.png)

``` r



## 設定を少し変えてもう一度実行する場合
sam_input <- res_sam$input # samに渡した引数のリスト
sam_input$rho.mode <- 1 # 引数の一部を変える（たとえば、Fの相関を1にする→今は選択率一定)
res_sam2 <- do.call(sam, sam_input) # do.callで再度計算
res_sam2$opt #収束している
#> $par
#>         logQ         logQ logSdLogFsta    logSdLogN  logSdLogObs  logSdLogObs 
#>     9.512568    -7.296031    -2.061695    -1.834227    -2.377935    -1.437568 
#>  logSdLogObs     rec_loga     rec_logb 
#>    -1.110561     2.382367   -19.094873 
#> 
#> $objective
#> [1] -19.7347
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 78
#> 
#> $evaluations
#> function gradient 
#>      105       79 
#> 
#> $message
#> [1] "relative convergence (4)"
res_sam2$rep #sdreportの結果
#> sdreport(.) result
#>                Estimate   Std. Error
#> logQ           9.512568 9.125255e-02
#> logQ          -7.296031 1.177387e-01
#> logSdLogFsta  -2.061695 2.787371e-01
#> logSdLogN     -1.834227 2.616238e-01
#> logSdLogObs   -2.377935 1.120733e-01
#> logSdLogObs   -1.437568 2.347150e-01
#> logSdLogObs   -1.110561 2.293440e-01
#> rec_loga       2.382367 5.989813e-02
#> rec_logb     -19.094873 2.338457e+04
#> Maximum gradient component: 4.651708e-06
res_sam2$rep$pdHess #Hessianが正定値を持つかどうか（持たない場合SEが計算できない）
#> [1] TRUE

## 2つのモデルのAICの比較
## - 選択率一定モデルのほうがAICが小さい（＝より予測力が高い）
c(res_sam$aic, res_sam2$aic)
#> [1]   9.390811 -21.469394

c(res_sam$rho, res_sam2$rho) #Fの相関係数
#> [1] 0 1

sam_input$SR <- "RW" #加入をRWに変更
res_sam3 <- do.call(sam,sam_input)
res_sam3$opt　
#> $par
#>         logQ         logQ logSdLogFsta    logSdLogN  logSdLogObs  logSdLogObs 
#>     9.505219    -7.303391    -2.047747    -1.254582    -2.374554    -1.444442 
#>  logSdLogObs 
#>    -1.108991 
#> 
#> $objective
#> [1] -14.68522
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 53
#> 
#> $evaluations
#> function gradient 
#>       84       53 
#> 
#> $message
#> [1] "relative convergence (4)"
res_sam3$rep
#> sdreport(.) result
#>               Estimate Std. Error
#> logQ          9.505219 0.09135781
#> logQ         -7.303391 0.11811402
#> logSdLogFsta -2.047747 0.28031293
#> logSdLogN    -1.254582 0.25257649
#> logSdLogObs  -2.374554 0.11286891
#> logSdLogObs  -1.444442 0.23289432
#> logSdLogObs  -1.108991 0.22963190
#> Maximum gradient component: 0.0001249072
res_sam3$rep$pdHess
#> [1] TRUE

c(res_sam$aic, res_sam2$aic, res_sam3$aic)
#> [1]   9.390811 -21.469394 -15.370439


## SAMだけの比較 (fishing_mortalityはFの平均？)
## - 信頼区間付きバージョン
## - VPAの結果も並列で示せるが、その場合vpa関数の実行にはTMB=TRUEオプションをつける必要あり
plot_samvpa(list(SAM=res_sam, SAM2=res_sam2,SAM3=res_sam3),
            what.plot = c("biomass","SSB","Recruitment","U","catch","fishing_mortality"))
```

![](sam_files/figure-html/comparison-3.png)

``` r


## F at age by year (境さんコード、関数化必要？時系列が長い場合への対応が必要)
res_samdata <- make_assess_result(res_sam)　# 境さん関数(make_assess_result)

(gg <- res_samdata %>% dplyr::filter(stat=="faa") %>%
    ggplot() +
    geom_ribbon(aes(x=Age, ymax=upper, ymin=lower, group=Model, fill=Model), alpha=0.5) +
    geom_line(aes(x=Age, y=Value, colour=Model, group=Model, size=Model), linetype=1) +
    facet_wrap(.~Year, scales = "free_y", nrow=10, ncol=2, dir="v") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_fill_manual("Model", values = c("#F8766D", "black")) +
    scale_colour_manual("Model", values = c("#F8766D", "black")) +
    scale_size_manual("Model", values = c(1, 0.5), guide = "none") +
    theme_bw() + 
    # labs(title="FAA", x=xlab, y=ylab) +
    theme(axis.text.x = element_text(size = 11, color = "black"),
          axis.text.y = element_text(size = 11, color = "black"),
          axis.line.x = element_line(linewidth = 0.3528), axis.line.y = element_line(linewidth = 0.3528),
          axis.minor.ticks.length = rel(0.5)))
```

![](sam_files/figure-html/comparison-4.png)

``` r


## F at age by year
(gg <- res_samdata %>% dplyr::filter(stat=="faa") %>%
    mutate(Age=factor(Age)) %>%
    ggplot() +
    geom_ribbon(aes(x=Year, ymax=upper, ymin=lower, group=Model, fill=Model), alpha=0.5) +
    geom_line(aes(x=Year, y=Value, colour=Model, group=Model), linetype=1) +
    facet_wrap(.~Age, scales = "free_y", nrow=10, ncol=2, dir="v") +
    scale_y_continuous(limits = c(0, NA)) +
    #  scale_fill_manual("Model", values = c("#F8766D", "black")) +
    #  scale_colour_manual("Model", values = c("#F8766D", "black")) +
    #  scale_size_manual("Model", values = c(1, 0.5), guide = "none") +
    theme_bw() + 
    # labs(title="FAA", x=xlab, y=ylab) +
    theme(axis.text.x = element_text(size = 11, color = "black"),
          axis.text.y = element_text(size = 11, color = "black"),
          axis.line.x = element_line(linewidth = 0.3528), axis.line.y = element_line(linewidth = 0.3528),
          axis.minor.ticks.length = rel(0.5)))
```

![](sam_files/figure-html/comparison-5.png)

## モデル選択

- rho.mode 0, 1, 2, 3のどれか？
- 再生産関係の関数と自己相関
- どの要素でsigmaを推定値し、どの年齢間でsigmaを共通にするか?
  - いろいろ試してAICの小さいものを選ぶ。（最新のcAICというものもあるようだが、未実装）

``` r

## varF(Fのプロセス誤差)とvarC（CAAの観測誤差）をどの年齢間で分けるかをstepAICで検討する

# いまは1歳魚以上のvarNを固定しているが、それも検討することも可能
# 検討する変数(VarF, VarC)を境界を入れる年齢Xについてのデータを生成(varとXを列名に使用)
# ここでは収束しやすさからRWの場合を使用する（BHを使うとlog(b)->-InfになりSEが発散する
grid = expand.grid(var=c("varF","varC"),X=1:6) %>%
  filter(!(var == "varF" & X == 6)) #5-6歳のFは同じと仮定しており6歳は不要なので除いておく

res_select <- select_sigma_grid(
  res_sam3, grid=grid, check_converge = TRUE, SEmax=Inf) # 

# stage 0 が元のモデル
# varがどの分散を分けたか
# whichはどの年齢で分けたか（1なら0歳と1歳魚以上で分ける）
# AIC最少の分け方をして、AICが下がらなくなるまで実行
# selectedはステージ内でAIC最少、bestは検討したすべてのモデルでAIC最少
res_select$tbl_sigma %>% knitr::kable()
```

| stage |       AIC | convergence | pdHess |     maxSE | var  | which | model    |
|------:|----------:|------------:|:-------|----------:|:-----|------:|:---------|
|     0 | -15.37044 |           0 | TRUE   | 0.2803129 | NA   |    NA | selected |
|     1 | -16.17782 |           0 | TRUE   | 0.3434833 | varF |     1 | best     |
|     1 |  70.30363 |           1 | FALSE  | 0.3066315 | varC |     1 | NA       |
|     1 | -15.20191 |           1 | TRUE   | 0.2926398 | varF |     2 | NA       |
|     1 | -13.64304 |           1 | TRUE   | 0.2808002 | varC |     2 | NA       |
|     1 | -13.38721 |           0 | TRUE   | 0.2896330 | varF |     3 | NA       |
|     1 | -13.37596 |           1 | TRUE   | 0.2804344 | varC |     3 | NA       |
|     1 | -14.33871 |           1 | TRUE   | 0.2845544 | varF |     4 | NA       |
|     1 | -13.40533 |           0 | TRUE   | 0.2806029 | varC |     4 | NA       |
|     1 | -15.81438 |           1 | TRUE   | 0.2892509 | varF |     5 | NA       |
|     1 | -13.66380 |           1 | TRUE   | 0.2789667 | varC |     5 | NA       |
|     1 | -13.78186 |           0 | TRUE   | 0.4950590 | varC |     6 | NA       |
|     2 |  68.47550 |           1 | FALSE  | 3.5031052 | varC |     1 | NA       |
|     2 | -14.57080 |           1 | TRUE   | 0.3559270 | varF |     2 | NA       |
|     2 | -15.13593 |           0 | TRUE   | 0.3199824 | varC |     2 | NA       |
|     2 | -14.64335 |           0 | TRUE   | 0.3457028 | varF |     3 | NA       |
|     2 | -14.25631 |           0 | TRUE   | 0.3436794 | varC |     3 | NA       |
|     2 | -14.26322 |           1 | TRUE   | 0.3400563 | varF |     4 | NA       |
|     2 | -14.25064 |           1 | TRUE   | 0.3432627 | varC |     4 | NA       |
|     2 | -15.45370 |           1 | TRUE   | 0.3382439 | varF |     5 | NA       |
|     2 | -14.71694 |           0 | TRUE   | 0.3587414 | varC |     5 | NA       |
|     2 | -15.15816 |           0 | TRUE   | 0.5952348 | varC |     6 | selected |

``` r



## Best modelの結果チェック
# FAAの誤差が0歳と1歳以上で分かれる
res_best <- res_select$bestres #AIC最少をベストモデルとする
cbind(
  "Age" = 0:6,
  "SD_caa" = res_best$sigma.logC,
  "SD_faa" = res_best$sigma.logF,
  "SD_naa" = res_best$sigma.logN
) %>% knitr::kable()
```

| Age |    SD_caa |    SD_faa |    SD_naa |
|----:|----------:|----------:|----------:|
|   0 | 0.0898453 | 0.0988359 | 0.2748462 |
|   1 | 0.0898453 | 0.1475825 | 0.0100000 |
|   2 | 0.0898453 | 0.1475825 | 0.0100000 |
|   3 | 0.0898453 | 0.1475825 | 0.0100000 |
|   4 | 0.0898453 | 0.1475825 | 0.0100000 |
|   5 | 0.0898453 | 0.1475825 | 0.0100000 |
|   6 | 0.0898453 | 0.1475825 | 0.0100000 |

## 再生産関係のプロット

``` r


## 単純なプロット
plot_SR_simple(res_sam2) #BHモデルを例に
```

![](sam_files/figure-html/plot_SR_simple-1.png)

``` r


## SAMの結果からfraysrで使える再生産関係のオブジェクトを作成するとfrasyrの関数が使える
SR_sam0 <- res_sam2 %>% make_SRres
biopar <- derive_biopar(res_sam2, derive_year=1998:2000)
## steepness, R0, B0, 決定論的なMSY管理基準値の計算
steepness_sam <- calc_steepness(SR="BH", rec_pars=SR_sam0$pars,M=biopar$M, waa=biopar$waa, maa=biopar$maa, faa=biopar$faa)
## 再生産関係に依存しない管理基準値の計算(frasyr::ref.F)
# FcurrentとPope=FALSEを明示的に引数に含める必要がある
Fcurrent <- rowMeans(res_sam2$faa[,as.character(1998:2000)])
refF_res <- ref.F(res_sam2,Fcurrent=Fcurrent,Pope=FALSE)
```

![](sam_files/figure-html/plot_SR_simple-2.png)

``` r

refF_res$summary
#>            Fcurrent Fmed Flow Fhigh      Fmax      F0.1 Fmean FpSPR.10.SPR
#> max       0.4565963   NA   NA    NA 0.5620096 0.3383640    NA    0.6752882
#> mean      0.4121442  NaN  NaN   NaN 0.5072950 0.3054225   NaN    0.6095454
#> Fref/Fcur 1.0000000   NA   NA    NA 1.2308677 0.7410572    NA    1.4789612
#>           FpSPR.20.SPR FpSPR.30.SPR FpSPR.40.SPR FpSPR.50.SPR FpSPR.60.SPR
#> max          0.4380875    0.3135295    0.2310909    0.1704472    0.1229955
#> mean         0.3954374    0.2830058    0.2085930    0.1538533    0.1110212
#> Fref/Fcur    0.9594637    0.6866668    0.5061164    0.3732996    0.2693746
#>           FpSPR.70.SPR FpSPR.80.SPR FpSPR.90.SPR
#> max         0.08434062   0.05193112   0.02417321
#> mean        0.07612962   0.04687535   0.02181983
#> Fref/Fcur   0.18471597   0.11373531   0.05294221
```

## モデル診断

### jitter analysis

- 初期値を変えて再推定し、目的関数（負の対数尤度）の値が変わらないかをチェックする

``` r

# ここでは10回だけ
jitterres = do_jitter(res_best,SD=0.1,nsim=10) #SDは初期値を乱数発生させる際のSD
#> 0: -16.089
#> 1: -16.089
#> 2: -16.089
#> 3: -16.089
#> 4: -16.089
#> 5: -16.089
#> 6: -16.089
#> 7: -16.089
#> 8: -16.089
#> 9: -16.089
#> 10: -16.089
knitr::kable(jitterres$resdat) #変わらない
```

|  ID | obj_value |
|----:|----------:|
|   0 | -16.08891 |
|   1 | -16.08891 |
|   2 | -16.08891 |
|   3 | -16.08891 |
|   4 | -16.08891 |
|   5 | -16.08891 |
|   6 | -16.08891 |
|   7 | -16.08891 |
|   8 | -16.08891 |
|   9 | -16.08891 |
|  10 | -16.08891 |

### 残差プロット

``` r


## 資源量指数 (frasyrのplot_residual_vpaと統合してもよい？）
## - これは通常のresidual
resid_sam <- index_plot(res_best); wrap_plots(resid_sam,nrow=1)
```

![](sam_files/figure-html/plot_residual-1.png)

``` r


## catch at age
caa_resid <- caa_plot(res_best); wrap_plots(caa_resid,ncol=2)
```

![](sam_files/figure-html/plot_residual-2.png)

``` r


## total catch weight の比較も必要かも
# 今関数はない

# catch at age (sakai ver.)

caa_obs0 <- dat$caa %>% rownames_to_column(var="Age") %>% pivot_longer(cols=-Age, names_to="Year", values_to="Value") %>% mutate(Data="Observation")
caa_est0 <- res_best$caa %>% as.data.frame %>% rownames_to_column(var="Age") %>% pivot_longer(cols=-Age, names_to="Year", values_to="Value") %>% mutate(Data="Estimation")
caa_obs_est <- bind_rows(caa_obs0, caa_est0)

xlab <- "年齢"
ylab <- "漁獲尾数（百万尾）"
scale <- 1000
g1 <- ggplot() +
  geom_bar(data=caa_obs0, stat="identity", aes(x=Age, y=Value/scale), col="black", fill="#FFFFB3") +
  geom_line(data=caa_est0, aes(x=Age, y=Value/scale, group=Data), col="#F8766D", linewidth=1.5) +
  facet_wrap(~Year, scales = "free_y", nrow=10, dir="v") +
  ylim(c(0, NA)) + scale_x_discrete(limit = c("0", "1", "2", "3", "4", "5", "6", "7", "8+")) +
  theme_bw() + labs(title="CAAプロット(棒グラフ：観測値, 折れ線：推定値)", x=xlab, y=ylab) +
  theme(axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.line.x = element_line(linewidth = 0.3528), axis.line.y = element_line(linewidth = 0.3528),
        axis.minor.ticks.length = rel(0.5), legend.position = "none")
```

### One-Step-Ahead (OSA) residuals

- SAMはランダム効果モデルなので、ハイパーパラメータとランダム効果の推定でデータを二重に使っていることになる
- ランダム効果でデータに対して過度に合わせることも可能であり（過剰適合）、通常の残差を診断に使うのは適切ではないと考えられている
- また、時系列解析なので残差は独立ではないので、独立とみなした診断（QQ
  plotなど）は不適
- あるサンプルを除いた時の予測値（One Step Ahead (OSA)
  prediction）とデータとの残差 (OSA residuals)
  を診断に使うのがより適切らしい
- ’TMB::oneStepPredict’を使って行うが、SAM用の関数’do_osa_resid()’を用意してある

``` r

osa.simple <- do_osa_resid(res_best)
#> [1] 90
#> [1] 89
#> [1] 88
#> [1] 87
#> [1] 86
#> [1] 85
#> [1] 84
#> [1] 83
#> [1] 82
#> [1] 81
#> [1] 80
#> [1] 79
#> [1] 78
#> [1] 77
#> [1] 76
#> [1] 75
#> [1] 74
#> [1] 73
#> [1] 72
#> [1] 71
#> [1] 70
#> [1] 69
#> [1] 68
#> [1] 67
#> [1] 66
#> [1] 65
#> [1] 64
#> [1] 63
#> [1] 62
#> [1] 61
#> [1] 60
#> [1] 59
#> [1] 58
#> [1] 57
#> [1] 56
#> [1] 55
#> [1] 54
#> [1] 53
#> [1] 52
#> [1] 51
#> [1] 50
#> [1] 49
#> [1] 48
#> [1] 47
#> [1] 46
#> [1] 45
#> [1] 44
#> [1] 43
#> [1] 42
#> [1] 41
#> [1] 40
#> [1] 39
#> [1] 38
#> [1] 37
#> [1] 36
#> [1] 35
#> [1] 34
#> [1] 33
#> [1] 32
#> [1] 31
#> [1] 30
#> [1] 29
#> [1] 28
#> [1] 27
#> [1] 26
#> [1] 25
#> [1] 24
#> [1] 23
#> [1] 22
#> [1] 21
#> [1] 20
#> [1] 19
#> [1] 18
#> [1] 17
#> [1] 16
#> [1] 15
#> [1] 14
#> [1] 13
#> [1] 12
#> [1] 11
#> [1] 10
#> [1] 9
#> [1] 8
#> [1] 7
#> [1] 6
#> [1] 5
#> [1] 4
#> [1] 3
#> [1] 2
#> [1] 1
gg_osa = plot_osa_resid(osa.simple); wrap_plots(gg_osa,ncol=2)
```

![](sam_files/figure-html/OSA%20residual-1.png)

### レトロスペクティブ解析

- ‘retro_sam()’ という関数で実行可能
- Retrospective forecastingも同時に実行可能でプロットもできる

``` r

retro_res = retro_sam(res_best,n=5) #nは年数

(g_retro = retro_plot(res_best,retro_res,start_year=1991,mohn_position="bottomleft"))
```

![](sam_files/figure-html/retro-1.png)

``` r


# retrospective forecastingもできる
(g_retro2 = retro_plot(res_best,retro_res,start_year=1991,mohn_position="bottomleft", forecast=TRUE))
```

![](sam_files/figure-html/retro-2.png)

### Hindcast cross validation

- レトロスペクティブ解析で最新のデータを順番に除いていき、除いたIndexを予測するHindcast
  cross validationを実行する
- Mean absolute scaled error (MASE) で予測精度を評価する

``` r

(g_hindcast <- plot_hindcastCV(res_best, retro_res, show_mase = TRUE, mase_position = "bottomleft", use_index = 1:2, log = FALSE))
```

![](sam_files/figure-html/hindcasting-1.png)

``` r


# MASEの結果を取り出す場合
res_mase = calc_mase(res_best, retro_res, log = FALSE)
res_mase$mase %>% knitr::kable()
```

| idx | index   | denominator |   numerator |      MASE |
|----:|:--------|------------:|------------:|----------:|
|   1 | Index 1 | 324.3231275 | 212.4539428 | 0.6550687 |
|   2 | Index 2 |   0.1126078 |   0.1329193 | 1.1803741 |

### Leave-one-out index analysis

- Indexを一つずつの除いて、推定値の頑健性、影響力のあるIndexを明らかにする（ジャックナイフ解析？）
- ’do_loo_index()’という関数を使って実行する

``` r

loo_res <- do_loo_index(res_best)
reslist <- list()
for ( j in 0:length(loo_res)) {
    if(j==0) {
      reslist[[j+1]] <- res_best
    } else {
      reslist[[j+1]] <- loo_res[[j]]
    }
}
(g_loo <- plot_samvpa(reslist,CI=0,
                    scenario_name = c("full",as.character(-1:-2))))
```

![](sam_files/figure-html/loo-index-1.png)

### プロファイル尤度

- Cachability
  qなどを変化させたときの対数尤度を調べて、収束しているか、どの程度尤度が変わるか（不確実性の程度、信頼区間）を調べる

``` r

# Catchabilityに対するプロファイル尤度
# 同じ名前のパラメータが複数ある場合は引数'which_param'で位置を指定できる
res_profile <- samprofile(res_best, "logQ", which_param=1, param_range=c(8, 11), length=25)

# 縦軸は負の対数尤度
res_profile$obj_tbl[-1,] %>%
  ggplot(aes(x=par, y=obj_value)) + geom_line(linewidth=0.5) +
  geom_point(data=res_profile$obj_tbl[1,],colour="red")
```

![](sam_files/figure-html/profile_likelihood-1.png)

``` r



# Mのプロファイル尤度 => 関数がないからこんな感じで手動で
Ms = seq(0.3,0.7,by=0.1)
scns = str_c("M:",Ms)
samres_Mlist = map(Ms, function(i) {
  res = res_best
  input = res$input
  input$dat$M[] <- i
  # p0_list = res$par_list; input$p0.list <- p0_list #初期値を利用するとき
  res.c = do.call(sam,input)
  return( res.c )
})

data.frame(
  M = Ms,
  obj_value = sapply(samres_Mlist, function(x) x$opt$objective)
) %>% ggplot(aes(x=M,y=obj_value)) + geom_line() + geom_point()
```

![](sam_files/figure-html/profile_likelihood-2.png)

``` r


plot_samvpa(samres_Mlist,scenario_name=scns,CI=0.)
```

![](sam_files/figure-html/profile_likelihood-3.png)

### ブートストラップ

- ’boo_sam()’で実行できる
- ノンパラメトリックブートストラップ（‘method=“n”’）は厳密ではないので、パラメトリックブートストラップを推奨（‘method=“p”’）
- delta法で求めたCIも図に加える場合は’draw_deltaCI=TRUE’

``` r


res_boot <- boo_sam(res_best, n=20,method="p",seed=1) #時間節約のため20回
#> -----1-----
#> -----2-----
#> -----3-----
#> -----4-----
#> -----5-----
#> -----6-----
#> -----7-----
#> -----8-----
#> -----9-----
#> -----10-----
#> -----11-----
#> -----12-----
#> -----13-----
#> -----14-----
#> -----15-----
#> -----16-----
#> -----17-----
#> -----18-----
#> -----19-----
#> -----20-----
(g1 = plot_boosam(res_best,res_boot, CI=0.95, draw_deltaCI = TRUE))
```

![](sam_files/figure-html/bootstrap-1.png)

### Self-test / cross test

- 真のモデルをVPA or SAMとして、疑似データを生成し、VPA or
  SAMのモデルを疑似データにフィットさせるというシミュレーションが可能
- 真のモデル (operating model, OM)
  と推定モデルが同じモデルであればselt-test, 異なればcross-test
- VPA or SAMのestimability (self-test), 異なる仮定に対する推定の頑健性
  (cross test) を調べられる

``` r

## 真のモデルはSAM
# pseudo dataを生成
pdata_sam <- popsim_vpasam(res_best, n=20)
# self-test (SAMで推定)
# 前述のブートストラップと同じ（はず）
fit_sam2sam <- fit2PSdata(res_best, PSdata=pdata_sam)
metric_sam2sam = sumup_popsim(res_best,fit_sam2sam) #真のモデルの結果を使うこと
# SSBの要約統計量を出力
metric_sam2sam$summary %>% filter(stat=="SSB") %>% knitr::kable()
```

| stat | year | age | RMSE | MAE | RMSRE | MARE | MedBias | MedRelBias | Median | CV | value_true | lower | upper |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SSB | 1991 | NA | 4.748965 | 3.472682 | 0.0395290 | 0.0289056 | 1.0141799 | 0.0084417 | 121.15305 | 0.0383032 | 120.13887 | 114.36772 | 130.68867 |
| SSB | 1992 | NA | 4.418594 | 3.281835 | 0.0366443 | 0.0272169 | 1.0258709 | 0.0085078 | 121.60662 | 0.0358219 | 120.58075 | 115.85956 | 130.86682 |
| SSB | 1993 | NA | 4.350415 | 3.373667 | 0.0359135 | 0.0278503 | -0.8631408 | -0.0071254 | 120.27265 | 0.0368647 | 121.13579 | 114.06842 | 130.44005 |
| SSB | 1994 | NA | 3.536073 | 2.867511 | 0.0317202 | 0.0257229 | -0.4448374 | -0.0039904 | 111.03228 | 0.0325566 | 111.47711 | 105.73816 | 117.71557 |
| SSB | 1995 | NA | 3.691662 | 2.958111 | 0.0358542 | 0.0287298 | 0.0720792 | 0.0007000 | 103.03514 | 0.0368093 | 102.96306 | 96.97200 | 110.11772 |
| SSB | 1996 | NA | 4.250619 | 3.570114 | 0.0479317 | 0.0402581 | -0.2056242 | -0.0023187 | 88.47510 | 0.0492152 | 88.68073 | 81.37054 | 96.19893 |
| SSB | 1997 | NA | 5.319414 | 4.260240 | 0.0725801 | 0.0581284 | -0.4107809 | -0.0056049 | 72.87944 | 0.0746003 | 73.29023 | 65.25896 | 83.31792 |
| SSB | 1998 | NA | 6.021707 | 4.745437 | 0.1031374 | 0.0812779 | -0.6461976 | -0.0110678 | 57.73910 | 0.1053884 | 58.38530 | 51.06980 | 70.83089 |
| SSB | 1999 | NA | 7.587664 | 5.989982 | 0.1377491 | 0.1087443 | -1.3693275 | -0.0248593 | 53.71388 | 0.1420179 | 55.08320 | 45.02114 | 70.36160 |
| SSB | 2000 | NA | 10.591320 | 8.474815 | 0.1947344 | 0.1558198 | -3.6841892 | -0.0677383 | 50.70436 | 0.2023717 | 54.38855 | 39.46957 | 76.36913 |

``` r

(g_sam2sam <- plot_popsim(res_best,fit_sam2sam))
```

![](sam_files/figure-html/cross_test-1.png)

``` r


# cross-test (VPAで推定)
# 前述のブートストラップと同じ（はず）
fit_vpa2sam <- fit2PSdata(res_vpa, PSdata=pdata_sam) #ここでは別のモデルを使う
metric_vpa2sam = sumup_popsim(res_best,fit_vpa2sam) #こっちでは真のモデルの結果を使うこと!
# SSBの要約統計量を出力
metric_vpa2sam$summary %>% filter(stat=="SSB") %>% knitr::kable()
```

| stat | year | age | RMSE | MAE | RMSRE | MARE | MedBias | MedRelBias | Median | CV | value_true | lower | upper |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SSB | 1991 | NA | 7.008667 | 5.639811 | 0.0583380 | 0.0469441 | 4.980693 | 0.0414578 | 125.11956 | 0.0363967 | 120.13887 | 119.12141 | 133.89929 |
| SSB | 1992 | NA | 6.803587 | 5.837697 | 0.0564235 | 0.0484132 | 5.359513 | 0.0444475 | 125.94026 | 0.0345226 | 120.58075 | 117.89615 | 133.15689 |
| SSB | 1993 | NA | 6.059764 | 5.088046 | 0.0500246 | 0.0420028 | 4.342049 | 0.0358445 | 125.47784 | 0.0374352 | 121.13579 | 116.63157 | 134.58281 |
| SSB | 1994 | NA | 5.614226 | 4.462884 | 0.0503621 | 0.0400341 | 3.639329 | 0.0326464 | 115.11644 | 0.0384839 | 111.47711 | 107.76149 | 122.85626 |
| SSB | 1995 | NA | 5.439524 | 4.135281 | 0.0528299 | 0.0401628 | 1.940047 | 0.0188422 | 104.90310 | 0.0395911 | 102.96306 | 100.25963 | 114.14214 |
| SSB | 1996 | NA | 6.117344 | 4.703847 | 0.0689817 | 0.0530425 | 2.256859 | 0.0254493 | 90.93758 | 0.0585842 | 88.68073 | 82.86073 | 100.84286 |
| SSB | 1997 | NA | 7.575826 | 5.668556 | 0.1033675 | 0.0773440 | 1.435129 | 0.0195815 | 74.72535 | 0.0950577 | 73.29023 | 65.45184 | 89.05097 |
| SSB | 1998 | NA | 8.263546 | 6.510862 | 0.1415347 | 0.1115154 | 2.368833 | 0.0405724 | 60.75413 | 0.1294029 | 58.38530 | 51.81063 | 77.45656 |
| SSB | 1999 | NA | 10.035738 | 8.013141 | 0.1821923 | 0.1454734 | 3.996335 | 0.0725509 | 59.07954 | 0.1646626 | 55.08320 | 46.85553 | 78.38818 |
| SSB | 2000 | NA | 12.793289 | 10.346382 | 0.2352202 | 0.1902309 | 4.315383 | 0.0793436 | 58.70394 | 0.2138959 | 54.38855 | 41.33169 | 83.14742 |

``` r

(g_vpa2sam <- plot_popsim(res_best,fit_vpa2sam))
```

![](sam_files/figure-html/cross_test-2.png)

``` r


# SAMの結果（res_best）の代わりに、VPAの結果をoperating model (真のモデル) として入れ替えて、VPAデータに対するself-test / cross-testもできます
```

## 将来予測

- 資源量推定の不確実性をどこまで考慮するか？
- 再生産関係は外部で推定しても良いのかも（難しいことはしなくてもよいようにもしておく）
