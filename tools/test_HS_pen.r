

setwd("~/git/frasam/tools")

library(devtools)

load_all("~/git/frasam")
library(frasyr)
library(TMB)

use_sam_tmb("HS_pen",overwrite=TRUE)

set.seed(123)
SSB = runif(50,100,1000)
R = exp(rnorm(50,log(SSB*0.1),0.5))

plot(R~SSB)

data_SR = list(R=R,SSB=SSB,years=1:50)

res = fit.SR(data_SR,SR="HS",AR=0,method="L2",length=50)
res$pars
range(SSB)

sum(res$resid^2)
max(SSB)

obj = MakeADFun(data=list(SSB=SSB,R=R),parameters=list(a=0.1,b=1000),random=NULL)
opt = nlminb(obj$par,obj$fn,obj$gr)
opt$par
(opt$par["b"]-max(SSB))^2
opt

res$pars

obj2 = MakeADFun(data=list(SSB=SSB,R=R),parameters=list(a=0.1,b=500),random=NULL)
opt2 = nlminb(obj2$par,obj2$fn,obj2$gr)
opt2$par
# (opt2$par["b"]-max(SSB))^2
opt$par
c(opt$objective,opt2$objective)


obj3 = MakeADFun(data=list(SSB=SSB,R=R),parameters=list(a=0.1,b=50),random=NULL)
opt3 = nlminb(obj3$par,obj3$fn,obj3$gr)
opt3$par
opt3
c(opt$objective,opt2$objective,opt3$objective)
cbind(opt$par,opt2$par,opt3$par)

sqrt(opt$objective/50)
res$pars


# obj4 = MakeADFun(data=list(SSB=SSB,R=R),parameters=list(a=0.1,b=1000),random=c("b"))
# obj4 = MakeADFun(data=list(SSB=SSB,R=R),parameters=list(a=0.1,b=500),random=c("b"))
obj4 = MakeADFun(data=list(SSB=SSB,R=R),parameters=list(a=0.1,b=50),random=c("b"))
opt4 = nlminb(obj4$par,obj4$fn,obj4$gr)
opt4
obj4$env$parList()["b"]
opt$par


obj5 = MakeADFun(data=list(SSB=SSB,R=R),parameters=list(a=0.1,b=500),random=c("a"))
opt5 = nlminb(obj5$par,obj5$fn,obj5$gr)
opt5
obj5$env$parList()[c("a","b")]
opt$par

