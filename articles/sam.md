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
library(frasam)
# devtools::load_all()
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
               est.method = "ml",
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
#> max       0.4565963   NA   NA    NA 0.5620096 0.3383543    NA    0.6752882
#> mean      0.4121442  NaN  NaN   NaN 0.5072950 0.3054137   NaN    0.6095454
#> Fref/Fcur 1.0000000   NA   NA    NA 1.2308677 0.7410360    NA    1.4789612
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

| stat | year | age |      RMSE |      MAE |     RMSRE |     MARE |    MedBias | MedRelBias |    Median |       CV | value_true |     lower |    upper |
|:-----|-----:|----:|----------:|---------:|----------:|---------:|-----------:|-----------:|----------:|---------:|-----------:|----------:|---------:|
| SSB  | 1991 |  NA |  96262.07 | 21528.51 |  801.2566 | 179.1968 |  1.6557709 |  0.0137821 | 121.79464 | 4.447057 |  120.13887 | 113.26746 | 226135.2 |
| SSB  | 1992 |  NA | 103773.23 | 23208.33 |  860.6119 | 192.4713 |  1.4633539 |  0.0121359 | 122.04410 | 4.448805 |  120.58075 | 112.58742 | 243771.1 |
| SSB  | 1993 |  NA | 110255.84 | 24658.07 |  910.1839 | 203.5572 |  0.5281819 |  0.0043602 | 121.66397 | 4.450090 |  121.13579 | 114.07817 | 258994.5 |
| SSB  | 1994 |  NA | 117114.74 | 26191.17 | 1050.5720 | 234.9466 | -0.5552580 | -0.0049809 | 110.92186 | 4.452982 |  111.47711 | 107.61641 | 275087.5 |
| SSB  | 1995 |  NA | 123716.25 | 27667.42 | 1201.5596 | 268.7121 | -0.7903540 | -0.0076761 | 102.17270 | 4.455436 |  102.96306 |  98.82600 | 290580.7 |
| SSB  | 1996 |  NA | 128514.97 | 28740.72 | 1449.1871 | 324.0920 | -1.2458636 | -0.0140489 |  87.43486 | 4.458343 |   88.68073 |  82.41123 | 301834.8 |
| SSB  | 1997 |  NA | 131688.45 | 29451.14 | 1796.8078 | 401.8427 | -1.4834311 | -0.0202405 |  71.80679 | 4.461028 |   73.29023 |  65.95243 | 309272.2 |
| SSB  | 1998 |  NA | 133517.90 | 29861.36 | 2286.8410 | 511.4533 | -2.5340113 | -0.0434015 |  55.85129 | 4.463367 |   58.38530 |  51.51799 | 313554.4 |
| SSB  | 1999 |  NA | 134641.20 | 30116.15 | 2444.3240 | 546.7392 | -5.3210389 | -0.0966000 |  49.76216 | 4.464039 |   55.08320 |  44.06890 | 316193.9 |
| SSB  | 2000 |  NA | 135429.14 | 30296.66 | 2490.0301 | 557.0411 | -7.8548661 | -0.1444213 |  46.53369 | 4.464263 |   54.38855 |  37.48370 | 318049.4 |

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

| stat | year | age |      RMSE |       MAE |     RMSRE |      MARE |    MedBias | MedRelBias |    Median |        CV | value_true |     lower |     upper |
|:-----|-----:|----:|----------:|----------:|----------:|----------:|-----------:|-----------:|----------:|----------:|-----------:|----------:|----------:|
| SSB  | 1991 |  NA |  6.605296 |  5.336627 | 0.0549805 | 0.0444205 |  4.3995653 |  0.0366207 | 124.53844 | 0.0383705 |  120.13887 | 117.28245 | 133.76663 |
| SSB  | 1992 |  NA |  6.836060 |  5.646956 | 0.0566928 | 0.0468313 |  4.1252020 |  0.0342111 | 124.70595 | 0.0386172 |  120.58075 | 117.24014 | 133.99612 |
| SSB  | 1993 |  NA |  7.200028 |  5.712401 | 0.0594377 | 0.0471570 |  4.2354347 |  0.0349644 | 125.37122 | 0.0393156 |  121.13579 | 119.24504 | 135.32616 |
| SSB  | 1994 |  NA |  6.118520 |  4.829909 | 0.0548859 | 0.0433265 |  4.5479614 |  0.0407973 | 116.02508 | 0.0344211 |  111.47711 | 110.90897 | 123.11240 |
| SSB  | 1995 |  NA |  5.392741 |  4.526218 | 0.0523755 | 0.0439596 |  4.3518814 |  0.0422664 | 107.31494 | 0.0345545 |  102.96306 | 100.93189 | 112.34283 |
| SSB  | 1996 |  NA |  5.019985 |  4.420800 | 0.0566074 | 0.0498507 |  3.5027798 |  0.0394988 |  92.18351 | 0.0470418 |   88.68073 |  83.87558 |  97.04532 |
| SSB  | 1997 |  NA |  5.269523 |  4.524783 | 0.0718994 | 0.0617379 |  2.4844352 |  0.0338986 |  75.77466 | 0.0670101 |   73.29023 |  67.34191 |  82.99795 |
| SSB  | 1998 |  NA |  6.380105 |  5.383204 | 0.1092759 | 0.0922014 |  0.2418703 |  0.0041427 |  58.62717 | 0.1070840 |   58.38530 |  51.25621 |  69.66393 |
| SSB  | 1999 |  NA | 10.802357 |  8.840768 | 0.1961098 | 0.1604984 | -0.1770578 | -0.0032144 |  54.90615 | 0.1921374 |   55.08320 |  44.84134 |  77.35118 |
| SSB  | 2000 |  NA | 17.519385 | 13.486490 | 0.3221153 | 0.2479656 | -0.9818225 | -0.0180520 |  53.40673 | 0.3075988 |   54.38855 |  39.36882 |  92.82609 |

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
