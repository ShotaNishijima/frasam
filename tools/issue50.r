
library(devtools)
library(frasyr)
library(tidyverse)

devtools::load_all()

use_sam_tmb(overwrite = TRUE)

data(sam_ex)

input <- sam_ex$input

input$p0.list <- sam_ex$par_list
input$min.age
input$max.age
input$abund
sam_ex$sigma.logN

input$varN.fix <- c(NA, 0.0001)
input$max.age[1:3] <- 1:3 # Nのage aggregate indexを作る

# load_all()

temp <- safe_call(sam, input)
res_index = index_plot(temp)

res_index$resid

naa <- temp$naa
pred.index <- temp$pred.index

temp$q[1]*(temp$naa[1,] + temp$naa[2,])
pred.index[1,]

resid = log(temp$input$dat$index[1,]) - log(temp$pred.index[1,])

as.numeric(resid) %>% mean

tmbdata = temp$data
obs = tmbdata$obs

obs_orig <- sam_ex$data$obs

obs_orig

obs_orig %>% as.data.frame() %>% filter(fleet > 1)

obs %>% as.data.frame() %>% filter(fleet > 1)

tmbdata$fleetTypes

sam_ex$data$fleetTypes

mean(log(temp$input$dat$index[1,]) - log(temp$pred.index[1,]), na.rm = T)

temp$rep

res_index = index_plot(temp)

res_index$abund
res_index$resid

res_index$index

temp


documdocument()
test()

install.packages("roxygen2")

data("dat_example")
dat_example <- dat
save(dat_example, file = "data/dat_example.rda")

plot_samvpa(sam_ex, CI=0.8)

plot_samvpa(list(res_rw,res_bh,res_ri), CI=0.8,
            scenario_name=rev(c("RW","BH","RI")))

sam_ex$caa %>% dimnames
res_rw$caa %>% dimnames

sam_ex$faa %>% dimnames
res_rw$faa %>% dimnames



# Kobe plot
# 管理基準値にラベルをつけるか？(0: つけない, 1:つける)
put_label_kobe <- 1
# Btarget, Blimit, Bbanのラベルを指定(kobe chart以外でも共通のラベルとなる）
label_name_kobe <- c("目標管理基準値案","限界管理基準値案","禁漁水準案")
# kobe plot を書く年の範囲 (0: 特に指定しない, 年数(1990:2000とか): その年のデータだけでKobe plotを書く
plot_year_kobe <- 1970:2023 #1989:2018
# 漁獲量曲線に過去の漁獲量を重ね書きする場合、過去の漁獲量をとる範囲 (0: 全年を指定)
past_year_range_yieldcurve <- 0 # c(1990:2000, 2005:2010) # 左のように、漁獲量を書く年をすべて指定する
# 再生産関係の点にラベルをつける年
SRplot_label_year <- c(1971,1978,1985,2004,2013,2018,2024)
# 再生産関係の図で示す予測区間の広さ(モデル平均の場合、現状では0.9で固定で調整できず）
predict_interval_SRplot <- 0.9
# HCRを重ね書きする場合のbeta (デフォルトは0.8, HCRを重ね書きしない場合は負の値を入れる)
# * HCRを入れる場合と入れない場合を2つ作る場合はベクトルで入力
# * ylabel_kobeとのすべての組み合わせのkobe chartが出力される
beta_kobe <- c(-1,0.8)
# HCRの図を書くときに設定するベータの値
beta_default <- 0.8

stopifnot(
  !is.null(sam_bh_prec$rep$jointPrecision),
  length(sam_bh_prec$obj$env$last.par.best) ==
    nrow(sam_bh_prec$rep$jointPrecision)
)

dir.create(
  "tests/testthat/testdata",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  sam_bh_prec,
  file = "tests/testthat/testdata/sam_bh_prec.rds",
  compress = "xz"
  )



