
library(devtools)
library(frasyr)
load_all()

res_sam2 <- readRDS("tools/res_sam2.RDS")
Fcurrent <- rowMeans(res_sam2$faa[,as.character(2020:2024)])

res_sam2$faa

refF_res <- ref.F(res_sam2,Fcurrent=Fcurrent,Pope=FALSE)

# 初期値を与えるとうまくいく
refF_res <- ref.F(res_sam2,Fcurrent=Fcurrent,Pope=FALSE, Fmax.init = 0.4)
# undebug(ref.F)

res_sam2$saa
refF_res$summary
refF_res$ypr.spr

