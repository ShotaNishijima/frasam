
library(tidyverse)
# library(frasam)

devtools::load_all("~/git/frasyr")
devtools::load_all()
use_sam_tmb(TmbFile = "sam2", overwrite = TRUE, compile = "always")

data("sam_ex")

sam_ex$input$p0.list <- sam_ex$par_list
res_biased = safe_call(sam, sam_ex$input)

sam_ex$caa -res_biased$caa #ほぼ同じ

input <- res_biased$input
# input$p0.list <- res_biased$par_list
input$bias.correct <- TRUE

res_corrected = safe_call(sam, input)

summary(res_biased$rep) %>% rownames %>% unique
summary(res_biased$rep) %>% rownames %>% unique

sdr_summary = summary(res_corrected$rep)

catch_sdr <- sdr_summary[rownames(sdr_summary) == "Catch_biomass",]

caa_est = res_corrected$caa
caa_est_biased = res_biased$caa
#

waa_sdr <- sdr_summary[rownames(sdr_summary) == "stockMeanWeight_true",]

waa_cpp <- matrix(waa_sdr[,1], ncol= ncol(caa_est), byrow=TRUE)

waa = input$dat$waa
waa - waa_cpp

sop_catch = colSums(caa_est*waa)
sop_catch_biased = colSums(caa_est_biased*waa)

catch_table = as.tibble(catch_sdr) %>%
  mutate(year = as.numeric(colnames(sam_ex$naa))) %>%
  select(year, everything()) %>%
  mutate(sop_catch = sop_catch,
         sop_catch_biased = sop_catch_biased)

catch_table %>% View

# 一致するようになった
