#################################################
#################################################

# code to run a4a on sardine NW of Iberian
# clean code for the selected runs

# 15/12/2015

#################################################
#################################################

#install.packages("ggplotFL", repos="http://flr-project.org/R")

# load libraries

library(FLa4a)
library(ggplotFL)
library(fields) # for the colour palette

# see functios on the library for help

library(help=FLa4a)

# set working directory

#wd <- "C:/IEO/wg's/WGHANSA/ISPRA_WK/data/a4a_NW/"
wd <- "C:/use/proyectos/a4a_sardine/a4a_nw_clean"
setwd(wd)

# load stock and indices object

#load("C:/IEO/wg's/WGHANSA/ISPRA_WK/data/NW.RData")
load("C:/use/proyectos/a4a_sardine/data/NW.RData")

# correct the NW.idx and save it

index(NW.idx[[1]])[,9]<-NA
index(NW.idx[[1]])[,17]<-NA
save(list=c('NW.idx','NW.stk'), file="C:/use/proyectos/a4a_sardine/data/NW_corrected.RData" )

# see range

range(NW.stk)

#################################################
#################################################

# understanding better the S4 classes

showClass("a4aFitSA")
showClass("SCAPars")
showClass("a4aStkParams")
showClass("submodel")

# get the ADMB template

getTPL(".")

#################################################
#################################################

# General names:

# S: south
# BB: Bay of Biscay
# NW: Northwest
# A: All
# IB: Iberian

# q1: smoother
# q2: constant
# q3: smoother with extra weight for survey
# q4: constant with extra weight for survey
# q5: constant but without the depm survey

# f: fit ll
# r: residuals ll
# s: stk object from ll
# mc: fit mcmc
# mcmc: fit mcmc as coda object
# smc: stk object from mcmc 

#################################################
#################################################

#===============================================

# q1: model with constant catchability across ages and years
# qmod <- list( ~ 1, ~ 1)
# fmod: bivariate tensor but adding smooth on age to decrease the number of df of the tensor on F

#===============================================


qmod <- list( ~ 1, ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
NW.q1f <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")
NW.q1r <- residuals(NW.q1f, NW.stk, NW.idx)
NW.q1s <- NW.stk + simulate(NW.q1f, 1000)
fitSumm(NW.q1f)
AIC(NW.q1f)
BIC(NW.q1f)

pdf("plots_NW_q1_ll.pdf", onefile=T)
wireframe(data ~ age + year, data = as.data.frame(harvest(NW.q1s)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(NW.q1f)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(NW.q1s)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(NW.q1s)), drape = TRUE, main="Catches")
plot(NW.q1r)
bubbles(NW.q1r)
plot(NW.q1f, NW.stk)
plot(NW.q1f, NW.idx[1])  
plot(NW.q1s)
dev.off()

NW.q1mc <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=500000,mcsave=250))
NW.q1mcmc <- as.mcmc(NW.q1mc)
NW.q1rmc <- residuals(NW.q1mc, NW.stk, NW.idx)
NW.q1smc <- NW.stk + NW.q1mc

plot(NW.q1smc)

#===============================================

# q2: model with smooth catchability across ages and years
# qmod <- list( ~ s(age, k=5), ~ 1)
# fmod: bivariate tensor but adding smooth on age to decrease the number of df of the tensor on F

#===============================================

# q2

qmod <- list( ~ s(age, k=5), ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
NW.q2f <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")
NW.q2r <- residuals(NW.q2f, NW.stk, NW.idx)
NW.q2s <- NW.stk + simulate(NW.q2f, 1000)
fitSumm(NW.q2f)
AIC(NW.q2f)
BIC(NW.q2f)

pdf("plots_NW_q2_ll.pdf", onefile=T)
wireframe(data ~ age + year, data = as.data.frame(harvest(NW.q2s)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(NW.q2f)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(NW.q2s)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(NW.q2s)), drape = TRUE, main="Catches")
plot(NW.q2r)
bubbles(NW.q2r)
plot(NW.q2f, NW.stk)
plot(NW.q2f, NW.idx[1])  
plot(NW.q2s)
dev.off()

NW.q2mc <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=500000,mcsave=250))
NW.q2mcmc <- as.mcmc(NW.q2mc)
NW.q2rmc <- residuals(NW.q2mc, NW.stk, NW.idx)
NW.q2smc <- NW.stk + NW.q2mc

plot(NW.q2smc)

#===============================================

# repeat q1 and q2 but giving more weight to surveys

#===============================================

# create a new index object with more weight to the acoustic survey

NW.idxw <- NW.idx
index.var(NW.idxw[[1]]) <- 0.5


#===============================================

# q3: model with constant catchability across ages and years
# but with more weights for acoustic survey
# qmod <- list( ~ 1, ~ 1)
# fmod: bivariate tensor but adding smooth on age to decrease the number of df of the tensor on F

#===============================================

# NW.q3

qmod <- list( ~ 1, ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
NW.q3f <- sca(NW.stk, NW.idxw, fmodel=fmod, qmodel=qmod, fit="assessment")
NW.q3r <- residuals(NW.q3f, NW.stk, NW.idxw)
NW.q3s <- NW.stk + simulate(NW.q3f, 1000)
fitSumm(NW.q3f)
AIC(NW.q3f)
BIC(NW.q3f)

pdf("plots_NW_q3_ll.pdf", onefile=T)
wireframe(data ~ age + year, data = as.data.frame(harvest(NW.q3s)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(NW.q3f)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(NW.q3s)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(NW.q3s)), drape = TRUE, main="Catches")
plot(NW.q3r)
bubbles(NW.q3r)
plot(NW.q3f, NW.stk)
plot(NW.q3f, NW.idxw[1])  
plot(NW.q3s)
dev.off()

NW.q3mc <- a4aSCA(NW.stk, NW.idxw, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=500000,mcsave=250))
NW.q3mcmc <- as.mcmc(NW.q3mc)
NW.q3rmc <- residuals(NW.q3mc, NW.stk, NW.idxw)
NW.q3smc <- NW.stk + NW.q3mc

plot(NW.q3smc)

#===============================================

# q4: model with smooth catchability across ages and years
# but with more weights for acoustic survey
# qmodel <- list( ~ s(age, k=5), ~ 1)
# fmod: bivariate tensor but adding smooth on age to decrease the number of df of the tensor on F

#===============================================

# NW.q4

qmod <- list( ~ s(age, k=5), ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
NW.q4f <- sca(NW.stk, NW.idxw, fmodel=fmod, qmodel=qmod, fit="assessment")
NW.q4r <- residuals(NW.q4f, NW.stk, NW.idxw)
NW.q4s <- NW.stk + simulate(NW.q4f, 1000)
fitSumm(NW.q4f)
AIC(NW.q4f)
BIC(NW.q4f)

pdf("plots_NW_q4_ll.pdf", onefile=T)
wireframe(data ~ age + year, data = as.data.frame(harvest(NW.q4s)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(NW.q4f)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(NW.q4s)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(NW.q4s)), drape = TRUE, main="Catches")
plot(NW.q4r)
bubbles(NW.q4r)
plot(NW.q4f, NW.stk)
plot(NW.q4f, NW.idxw[1])  
plot(NW.q4s)
dev.off()

NW.q4mc <- a4aSCA(NW.stk, NW.idxw, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=500000,mcsave=250))
NW.q4mcmc <- as.mcmc(NW.q4mc)
NW.q4rmc <- residuals(NW.q4mc, NW.stk, NW.idxw)
NW.q4smc <- NW.stk + NW.q4mc

plot(NW.q4smc)

#===============================================

# q5: model with smooth catchability across ages and years
# but without depm, using only acoustics
# qmod <- list( ~ s(age, k=5), ~ 1)
# fmod: bivariate tensor but adding smooth on age to decrease the number of df of the tensor on F

#===============================================

# NW.q5

qmod <- list( ~ s(age, k=5), ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
NW.q5f <- sca(NW.stk, NW.idx[1], fmodel=fmod, qmodel=qmod, fit="assessment")
NW.q5r <- residuals(NW.q5f, NW.stk, NW.idx[1])
NW.q5s <- NW.stk + simulate(NW.q5f, 1000)
fitSumm(NW.q5f)
AIC(NW.q5f)
BIC(NW.q5f)

pdf("plots_NW_q5_ll.pdf", onefile=T)
wireframe(data ~ age + year, data = as.data.frame(harvest(NW.q5s)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(NW.q5f)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(NW.q5s)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(NW.q5s)), drape = TRUE, main="Catches")
plot(NW.q5r)
bubbles(NW.q5r)
plot(NW.q5f, NW.stk)
plot(NW.q5f, NW.idx[1])  
plot(NW.q5s)
dev.off()

NW.q5mc <- a4aSCA(NW.stk, NW.idx[1], fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=500000,mcsave=250))
NW.q5mcmc <- as.mcmc(NW.q5mc)
NW.q5rmc <- residuals(NW.q5mc, NW.stk, NW.idx[1])
NW.q5smc <- NW.stk + NW.q5mc

plot(NW.q5smc)

#################################################
#################################################

# comparison of gcv, AIC and BIC

cbind(gcv=c(fitSumm(NW.q1f)['gcv',], fitSumm(NW.q2f)['gcv',], fitSumm(NW.q3f)['gcv',], fitSumm(NW.q4f)['gcv',], fitSumm(NW.q5f)['gcv',]),
      aic=AIC(NW.q1f, NW.q2f, NW.q3f, NW.q4f, NW.q5f),
      bic=BIC(NW.q1f, NW.q2f, NW.q3f, NW.q4f, NW.q5f))

#################################################
#################################################

# save the objects for the general analysis

save(list=c("NW.stk","NW.idx", "NW.idxw", paste("NW.", rep(c("q1","q2","q3","q4","q5"), each=6), rep(c("f","r","s","mc","mcmc","smc"),5), sep="")),
     file="NWresults.RData")

#################################################
#################################################

# comparison of loglik and MCMC intervals for all runs:

pdf("comparison_ll_mcmc.pdf")
plot(FLStocks(q1.ll=NW.q1s, q1.mc=NW.q1smc))
plot(FLStocks(q2.ll=NW.q2s, q2.mc=NW.q2smc))
plot(FLStocks(q3.ll=NW.q3s, q3.mc=NW.q3smc))
plot(FLStocks(q4.ll=NW.q4s, q4.mc=NW.q4smc))
plot(FLStocks(q5.ll=NW.q5s, q5.mc=NW.q5smc))
dev.off()

# comparison between cases using ll

pdf("comparison_ll.pdf")
plot(FLStocks(q1.ll=NW.q1s, q2.ll=NW.q2s, q3.ll=NW.q3s, q4.ll=NW.q4s, q5.ll=NW.q5s))
dev.off()

# comparison between cases using mcmc

pdf("comparison_mcmc.pdf")
plot(FLStocks(q1.mc=NW.q1smc, q2.mc=NW.q2smc, q3.mc=NW.q3smc, q4.mc=NW.q4smc, q5.mc=NW.q5smc))
dev.off()

#################################################
#################################################

# save the session

save.image("main.RData")

#################################################
#################################################

