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
#>    -1.082575     2.365818   -21.005623 
#> 
#> $objective
#> [1] -4.304594
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 1
#> 
#> $evaluations
#> function gradient 
#>        2        1 
#> 
#> $message
#> [1] "relative convergence (4)"

## 推定結果の確認: rep=推定パラメータとSE
## - 推定値(Estimate)に対してStd.Errorが異常に大きくないか
## - Maximum gradient componentが十分小さい値か（1e-3以下なら十分OK、多数のパラメータを推定する難しいモデルなら上限0.1以下くらいでも？)
res_sam$rep
#> sdreport(.) result
#>                Estimate   Std. Error
#> logQ           9.560474 8.811787e-02
#> logQ          -7.232696 1.118901e-01
#> logSdLogFsta  -1.877590 1.476097e-01
#> logSdLogN     -1.761564 2.699046e-01
#> logSdLogObs   -3.101510 7.511985e-01
#> logSdLogObs   -1.367962 2.342732e-01
#> logSdLogObs   -1.082575 2.288350e-01
#> rec_loga       2.365818 6.530280e-02
#> rec_logb     -21.005623 4.833104e+04
#> Maximum gradient component: 8.678031e-06

## 推定結果の確認: AIC（単独ではあまり意味がない。モデルを比較するときに）
res_sam$aic
#> [1] 9.390811

## 推定結果の確認: 分散パラメータのsigma。ノーマルスケールで、大きさを確認。
## - すごく小さい値で推定されている場合には推定の意味がないので推定しないなどする
## - どのsigmaを共通にするか？を決める
res_sam[c("sigma","sigma.logC","sigma.logFsta","sigma.logN")]
#> $sigma
#> [1] 0.2546254 0.3387221
#> 
#> $sigma.logC
#> [1] 0.04498125 0.04498125 0.04498125 0.04498125 0.04498125 0.04498125 0.04498125
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
| logQ         |   9.560474 | 8.811790e-02 |  8.0e-06 |   1.419257e+04 |
| logQ         |  -7.232697 | 1.118901e-01 |  7.9e-06 |   7.226000e-04 |
| logSdLogFsta |  -1.877590 | 1.476097e-01 |  8.2e-06 |   1.529583e-01 |
| logSdLogN    |  -1.761564 | 2.699046e-01 | -1.0e-07 |   1.717760e-01 |
| logSdLogObs  |  -3.101510 | 7.511985e-01 |  1.4e-06 |   4.498120e-02 |
| logSdLogObs  |  -1.367962 | 2.342732e-01 |  3.1e-06 |   2.546254e-01 |
| logSdLogObs  |  -1.082575 | 2.288350e-01 | -1.0e-07 |   3.387221e-01 |
| rec_loga     |   2.365818 | 6.530280e-02 |  8.7e-06 |   1.065275e+01 |
| rec_logb     | -21.005623 | 4.833104e+04 |  0.0e+00 |   0.000000e+00 |

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
#>     9.512569    -7.296031    -2.061696    -1.834225    -2.377936    -1.437569 
#>  logSdLogObs     rec_loga     rec_logb 
#>    -1.110560     2.382367   -19.991679 
#> 
#> $objective
#> [1] -19.7347
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 81
#> 
#> $evaluations
#> function gradient 
#>      106       82 
#> 
#> $message
#> [1] "relative convergence (4)"
res_sam2$rep #sdreportの結果
#> sdreport(.) result
#>                Estimate   Std. Error
#> logQ           9.512569 9.125957e-02
#> logQ          -7.296031 1.177465e-01
#> logSdLogFsta  -2.061696 2.792095e-01
#> logSdLogN     -1.834225 2.616243e-01
#> logSdLogObs   -2.377936 1.120733e-01
#> logSdLogObs   -1.437569 2.347180e-01
#> logSdLogObs   -1.110560 2.293442e-01
#> rec_loga       2.382367 5.989870e-02
#> rec_logb     -19.991679 3.667806e+04
#> Maximum gradient component: 4.786195e-05
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
#>     9.505217    -7.303391    -2.047742    -1.254581    -2.374555    -1.444445 
#>  logSdLogObs 
#>    -1.108989 
#> 
#> $objective
#> [1] -14.68522
#> 
#> $convergence
#> [1] 1
#> 
#> $iterations
#> [1] 4
#> 
#> $evaluations
#> function gradient 
#>       15        4 
#> 
#> $message
#> [1] "false convergence (8)"
res_sam3$rep
#> sdreport(.) result
#>               Estimate Std. Error
#> logQ          9.505217 0.09135759
#> logQ         -7.303391 0.11811426
#> logSdLogFsta -2.047742 0.28031350
#> logSdLogN    -1.254581 0.25257667
#> logSdLogObs  -2.374555 0.11286874
#> logSdLogObs  -1.444445 0.23289363
#> logSdLogObs  -1.108989 0.22963229
#> Maximum gradient component: 0.0001080425
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

| stage |       AIC | convergence | pdHess |     maxSE | var  | which | model |
|------:|----------:|------------:|:-------|----------:|:-----|------:|:------|
|     0 | -15.37044 |           1 | TRUE   | 0.2803135 | NA   |    NA | NA    |
|     1 | -16.17782 |           1 | TRUE   | 0.3466419 | varF |     1 | NA    |
|     1 |  70.29138 |           1 | FALSE  | 0.5953598 | varC |     1 | NA    |
|     1 | -15.20191 |           1 | TRUE   | 0.2895434 | varF |     2 | NA    |
|     1 | -13.64304 |           0 | TRUE   | 0.2807999 | varC |     2 | NA    |
|     1 | -13.38721 |           1 | TRUE   | 0.2894631 | varF |     3 | NA    |
|     1 | -13.37596 |           1 | TRUE   | 0.2804328 | varC |     3 | NA    |
|     1 | -14.33871 |           1 | TRUE   | 0.2848443 | varF |     4 | NA    |
|     1 | -13.40533 |           0 | TRUE   | 0.2806029 | varC |     4 | NA    |
|     1 | -15.81438 |           1 | TRUE   | 0.2884235 | varF |     5 | NA    |
|     1 | -13.66380 |           0 | TRUE   | 0.2789662 | varC |     5 | NA    |
|     1 | -13.78186 |           0 | TRUE   | 0.4950591 | varC |     6 | best  |

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
|   0 | 0.0930559 | 0.1290259 | 0.2851952 |
|   1 | 0.0930559 | 0.1290259 | 0.0100000 |
|   2 | 0.0930559 | 0.1290259 | 0.0100000 |
|   3 | 0.0930559 | 0.1290259 | 0.0100000 |
|   4 | 0.0930559 | 0.1290259 | 0.0100000 |
|   5 | 0.0930559 | 0.1290259 | 0.0100000 |
|   6 | 0.0930559 | 0.1290259 | 0.0100000 |

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
#> max       0.4565964   NA   NA    NA 0.5620096 0.3383542    NA    0.6752882
#> mean      0.4121444  NaN  NaN   NaN 0.5072950 0.3054137   NaN    0.6095454
#> Fref/Fcur 1.0000000   NA   NA    NA 1.2308674 0.7410357    NA    1.4789608
#>           FpSPR.20.SPR FpSPR.30.SPR FpSPR.40.SPR FpSPR.50.SPR FpSPR.60.SPR
#> max          0.4380875    0.3135295    0.2310909    0.1704472    0.1229955
#> mean         0.3954374    0.2830058    0.2085930    0.1538533    0.1110212
#> Fref/Fcur    0.9594635    0.6866666    0.5061163    0.3732995    0.2693746
#>           FpSPR.70.SPR FpSPR.80.SPR FpSPR.90.SPR
#> max         0.08434062   0.05193112   0.02417321
#> mean        0.07612962   0.04687536   0.02181983
#> Fref/Fcur   0.18471592   0.11373528   0.05294219
```

## モデル診断

### jitter analysis

- 初期値を変えて再推定し、目的関数（負の対数尤度）の値が変わらないかをチェックする

``` r

# ここでは10回だけ
jitterres = do_jitter(res_best,SD=0.1,nsim=10) #SDは初期値を乱数発生させる際のSD
#> 0: -14.685
#> 1: -14.685
#> 2: -14.685
#> 3: -14.685
#> 4: -14.685
#> 5: -14.685
#> 6: -14.685
#> 7: -14.685
#> 8: -14.685
#> 9: -14.685
#> 10: -14.685
knitr::kable(jitterres$resdat) #変わらない
```

|  ID | obj_value |
|----:|----------:|
|   0 | -14.68522 |
|   1 | -14.68522 |
|   2 | -14.68522 |
|   3 | -14.68522 |
|   4 | -14.68522 |
|   5 | -14.68522 |
|   6 | -14.68522 |
|   7 | -14.68522 |
|   8 | -14.68522 |
|   9 | -14.68522 |
|  10 | -14.68522 |

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
|   1 | Index 1 | 324.3231275 | 217.3824785 | 0.6702651 |
|   2 | Index 2 |   0.1126078 |   0.1314383 | 1.1672222 |

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
- Ridge VPAを推定モデルとして使用する場合、`reest_lambda = TRUE` とし、
  `target_retro`（例：`"SSB"`, `"B"`）を指定すると、各疑似データセットで
  レトロスペクティブバイアスに基づいてλを選択してから推定する。
  デフォルトの `FALSE` では、元の推定結果のλを固定して使用する。
  ηが設定されている場合のλ再推定は未対応。詳細は
  [`?fit2PSdata`](https://shotanishijima.github.io/frasam/reference/fit2PSdata.md)
  を参照。

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

| stat | year | age | RMSE | MAE | RMSRE | MedAbsRelBias | MARE | MedBias | MedRelBias | Median | CV | value_true | lower | upper |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SSB | 1991 | NA | 4.760883 | 3.414062 | 0.0399543 | 0.0229347 | 0.0286515 | 0.8649042 | 0.0072584 | 120.02326 | 0.0385574 | 119.15836 | 113.71365 | 129.88430 |
| SSB | 1992 | NA | 4.322414 | 3.067661 | 0.0362373 | 0.0190846 | 0.0257180 | 0.7651554 | 0.0064147 | 120.04603 | 0.0352502 | 119.28088 | 114.90947 | 129.71126 |
| SSB | 1993 | NA | 4.247140 | 3.274844 | 0.0355473 | 0.0187200 | 0.0274095 | -0.4367555 | -0.0036555 | 119.04187 | 0.0364825 | 119.47863 | 112.73112 | 128.82634 |
| SSB | 1994 | NA | 3.382955 | 2.708013 | 0.0307099 | 0.0162689 | 0.0245829 | -0.4867405 | -0.0044185 | 109.67168 | 0.0315132 | 110.15842 | 104.48736 | 116.39839 |
| SSB | 1995 | NA | 3.577191 | 2.817012 | 0.0351867 | 0.0217436 | 0.0277093 | -0.0277509 | -0.0002730 | 101.63538 | 0.0361111 | 101.66313 | 95.47938 | 108.51959 |
| SSB | 1996 | NA | 4.095084 | 3.338722 | 0.0470797 | 0.0309660 | 0.0383841 | 0.6931172 | 0.0079685 | 87.67508 | 0.0481243 | 86.98196 | 80.09680 | 94.42109 |
| SSB | 1997 | NA | 4.944481 | 3.746109 | 0.0695939 | 0.0445430 | 0.0527267 | 0.4225220 | 0.0059470 | 71.47015 | 0.0710478 | 71.04763 | 63.69307 | 80.96735 |
| SSB | 1998 | NA | 5.490464 | 4.028400 | 0.0992943 | 0.0724257 | 0.0728531 | -0.1234117 | -0.0022319 | 55.17145 | 0.0994167 | 55.29487 | 49.38485 | 67.22036 |
| SSB | 1999 | NA | 6.789493 | 4.955005 | 0.1357858 | 0.0800636 | 0.0990971 | 0.4117971 | 0.0082357 | 50.41330 | 0.1371719 | 50.00150 | 42.56497 | 65.50850 |
| SSB | 2000 | NA | 9.559206 | 6.813286 | 0.2043509 | 0.1105765 | 0.1456503 | 0.8893981 | 0.0190130 | 47.66780 | 0.2070169 | 46.77840 | 35.70513 | 68.96615 |

``` r

(g_sam2sam <- plot_popsim(res_best,fit_sam2sam))
```

![](sam_files/figure-html/cross_test-1.png)

``` r


# cross-test (VPAで推定)
# 前述のブートストラップと同じ（はず）
fit_vpa2sam <- fit2PSdata(res_vpa, PSdata=pdata_sam, reest_lambda = FALSE) #ここでは別のモデルを使う
metric_vpa2sam = sumup_popsim(res_best,fit_vpa2sam) #こっちでは真のモデルの結果を使うこと!
# SSBの要約統計量を出力
metric_vpa2sam$summary %>% filter(stat=="SSB") %>% knitr::kable()
```

| stat | year | age | RMSE | MAE | RMSRE | MedAbsRelBias | MARE | MedBias | MedRelBias | Median | CV | value_true | lower | upper |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SSB | 1991 | NA | 7.186206 | 5.799497 | 0.0603080 | 0.0434618 | 0.0486705 | 5.178836 | 0.0434618 | 124.33720 | 0.0372999 | 119.15836 | 118.21228 | 133.26274 |
| SSB | 1992 | NA | 6.984090 | 6.009250 | 0.0585516 | 0.0469380 | 0.0503790 | 5.598801 | 0.0469380 | 124.87968 | 0.0353237 | 119.28088 | 116.65272 | 132.18572 |
| SSB | 1993 | NA | 6.146328 | 5.178209 | 0.0514429 | 0.0385754 | 0.0433400 | 4.442901 | 0.0371857 | 123.92153 | 0.0378161 | 119.47863 | 115.11378 | 133.05367 |
| SSB | 1994 | NA | 5.669393 | 4.554982 | 0.0514658 | 0.0384308 | 0.0413494 | 3.900287 | 0.0354062 | 114.05870 | 0.0381415 | 110.15842 | 106.65218 | 121.56909 |
| SSB | 1995 | NA | 5.484904 | 4.246310 | 0.0539517 | 0.0311311 | 0.0417684 | 2.454582 | 0.0241443 | 104.11772 | 0.0383099 | 101.66313 | 99.44525 | 112.68431 |
| SSB | 1996 | NA | 6.053328 | 4.704203 | 0.0695929 | 0.0451323 | 0.0540825 | 2.946194 | 0.0338713 | 89.92815 | 0.0560435 | 86.98196 | 81.85830 | 98.93936 |
| SSB | 1997 | NA | 7.351547 | 5.431209 | 0.1034735 | 0.0438145 | 0.0764446 | 2.147906 | 0.0302319 | 73.19554 | 0.0910254 | 71.04763 | 64.21303 | 86.33457 |
| SSB | 1998 | NA | 7.898897 | 6.057660 | 0.1428505 | 0.0853502 | 0.1095519 | 3.411076 | 0.0616888 | 58.70594 | 0.1229980 | 55.29487 | 50.17621 | 73.63776 |
| SSB | 1999 | NA | 9.476641 | 7.419899 | 0.1895271 | 0.1206661 | 0.1483935 | 4.950187 | 0.0990008 | 54.95169 | 0.1585826 | 50.00150 | 43.96213 | 72.16323 |
| SSB | 2000 | NA | 12.038027 | 9.527986 | 0.2573416 | 0.1802013 | 0.2036835 | 5.661384 | 0.1210256 | 52.43978 | 0.2148972 | 46.77840 | 36.65093 | 74.20437 |

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
