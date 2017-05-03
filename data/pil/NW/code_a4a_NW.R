#################################################
#################################################

# code to run a4a on sardine NW of Iberian

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
wd <- "C:/use/proyectos/a4a_sardine/a4a_nw"
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

#===============================================

# models with constant catchability across ages and years
# qmodel <- list( ~ 1, ~ 1)

#===============================================

# fit1
# default options for fmodel, srmodel, n1model and vmodel

qmod <- list( ~ 1, ~ 1)
fit1c <- sca(NW.stk, NW.idx, qmodel=qmod, fit="assessment")

fit <- fit1c

pdf("plots_fit1c.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#-----------------------------------------------

# fit2c
# increase the number of knots for age and year in F

qmod <- list( ~ 1, ~ 1)
fmod <- ~ te(age, year, k=c(6,15))
fit2c <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit2c
pdf("plots_fit2c.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#-----------------------------------------------

# fit3c
# intermediate number of knots for age and year in F

qmod <- list( ~ 1, ~ 1)
fmod <- ~ te(age, year, k=c(5,12))
fit3c <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit3c
pdf("plots_fit3c.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#-----------------------------------------------

# fit4c
# separable F with smooths

qmod <- list( ~ 1, ~ 1)
fmod <- ~ s(age, k=5) + s(year, k=12)
fit4c <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit4c
pdf("plots_fit4c.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#-----------------------------------------------

# fit5c
# separable F with factors doesn't converge,so use separable factor and smooth

qmod <- list( ~ 1, ~ 1)
fmod <- ~ factor(age) + s(year, k=12)  
fit5c <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit5c
pdf("plots_fit5c.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#-----------------------------------------------

# fit6c
# to decrease the number of df of the tenson on F

qmod <- list( ~ 1, ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6c <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit6c
pdf("plots_fit6c.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#===============================================

# models with smooth catchability across ages and years
# qmodel <- list( ~ s(age, k=5), ~ 1)

#===============================================

# fit1s
# default options for fmodel, srmodel, n1model and vmodel

qmod <- list( ~ s(age, k=5), ~ 1)
fit1s <- sca(NW.stk, NW.idx, qmodel=qmod, fit="assessment")

fit <- fit1s

pdf("plots_fit1s.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#-----------------------------------------------

# fit2s
# increase the number of knots for age and year in F

qmod <- list( ~ s(age, k=5), ~ 1)
fmod <- ~ te(age, year, k=c(6,14))
fit2s <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit2s
pdf("plots_fit2s.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#-----------------------------------------------

# fit3s
# intermediate number of knots for age and year in F

qmod <- list( ~ s(age, k=5), ~ 1)
fmod <- ~ te(age, year, k=c(5,12))
fit3s <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit3s
pdf("plots_fit3s.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#-----------------------------------------------

# fit4s
# separable F with smooths

qmod <- list( ~ s(age, k=5), ~ 1)
fmod <- ~ s(age, k=5) + s(year, k=12)
fit4s <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit4s
pdf("plots_fit4s.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#-----------------------------------------------

# fit6s
# to decrease the number of df of the tensor on F

qmod <- list( ~ s(age, k=5), ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6s <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit6s
pdf("plots_fit6s.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#===============================================

# models with catchability for age 1 and age 2+
# qmodel <- list( ~ factor(replace(age, age >1, 2)), ~ 1)

#===============================================

qmod <- list( ~ factor(replace(age, age >1, 2)), ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6q <- sca(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit6q
pdf("plots_fit6q.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#===============================================

# repeat fit6c, fit6s and fit6q but giving more weight to surveys

#===============================================

# create a new index object with more weight to the acoustic survey

NW.idxw <- NW.idx
index.var(NW.idxw[[1]]) <- 0.5

# fit6cw

qmod <- list( ~ 1, ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6cw <- sca(NW.stk, NW.idxw, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit6cw
pdf("plots_fit6cw.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idxw)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idxw[1])  
plot(stk)
dev.off()

# fit6sw

qmod <- list( ~ s(age, k=5), ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6sw <- sca(NW.stk, NW.idxw, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit6sw
pdf("plots_fit6sw.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idxw)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idxw[1])  
plot(stk)
dev.off()

# fit6qw

qmod <- list( ~ factor(replace(age, age >1, 2)), ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6qw <- sca(NW.stk, NW.idxw, fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit6qw
pdf("plots_fit6qw.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idxw)
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idxw[1])  
plot(stk)
dev.off()

#===============================================

# repeat fit6c but using only the depm survey

#===============================================

# fit6co

qmod <- list( ~ 1, ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6co <- sca(NW.stk, NW.idx[1], fmodel=fmod, qmodel=qmod, fit="assessment")

fit <- fit6co
pdf("plots_fit6co.pdf", onefile=T)
res <- residuals(fit, NW.stk, NW.idx[1])
stk <- NW.stk + fit
fitSumm(fit)
AIC(fit)
BIC(fit)
wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality") 
wireframe(data ~ age + year, data = as.data.frame(predict(fit)$qmodel[[1]]), drape = TRUE, screen = list(x = -90, y=-45), main="Catchability")
wireframe(data ~ age + year, data = as.data.frame(stock.n(stk)), drape = TRUE, main="Population", screen = list(x = -90, y=-45))
wireframe(data ~ age + year, data = as.data.frame(catch.n(stk)), drape = TRUE, main="Catches")
plot(res)
bubbles(res)
plot(fit, NW.stk)
plot(fit, NW.idx[1])  
plot(stk)
dev.off()

#################################################
#################################################

# comparison of gcv, AIC and BIC

cbind(gcv=c(fitSumm(fit6c)['gcv',], fitSumm(fit6s)['gcv',], fitSumm(fit6q)['gcv',], fitSumm(fit6cw)['gcv',], fitSumm(fit6sw)['gcv',], fitSumm(fit6qw)['gcv',]),
      aic=AIC(fit6c, fit6s, fit6q, fit6cw, fit6sw, fit6qw),
      bic=BIC(fit6c, fit6s, fit6q, fit6cw, fit6sw, fit6qw))

#################################################
#################################################

# compute confidence intervals using the hessian

stk6c.ll <- NW.stk + simulate(fit6c, 1000)
stk6s.ll <- NW.stk + simulate(fit6s, 1000)
stk6q.ll <- NW.stk + simulate(fit6q, 1000)
stk6cw.ll <- NW.stk + simulate(fit6cw, 1000)
stk6sw.ll <- NW.stk + simulate(fit6sw, 1000)
stk6qw.ll <- NW.stk + simulate(fit6qw, 1000)

#################################################
#################################################

# compute confidence intervals using mcmcm

# fit6c

qmod <- list( ~ 1, ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6c.mc <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=20000,mcsave=20))

# fit6s

qmod <- list( ~ s(age, k=5), ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6s.mc <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=20000,mcsave=20))

# fit6q

qmod <- list( ~ factor(replace(age, age >1, 2)), ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6q.mc <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=20000,mcsave=20))

# fit6cw

qmod <- list( ~ 1, ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6cw.mc <- a4aSCA(NW.stk, NW.idxw, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=20000,mcsave=20))

# fit6sw

qmod <- list( ~ s(age, k=5), ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6sw.mc <- a4aSCA(NW.stk, NW.idxw, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=20000,mcsave=20))

# fit6qw

qmod <- list( ~ factor(replace(age, age >1, 2)), ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6qw.mc <- a4aSCA(NW.stk, NW.idxw, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=20000,mcsave=20))

# fit6co.mc

qmod <- list( ~ 1, ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
fit6co.mc <- a4aSCA(NW.stk, NW.idx[1], fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=20000,mcsave=20))

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

NW.q1f <- fit6s
NW.q1r <- residuals(NW.q1f, NW.stk, NW.idx)
NW.q1s <- NW.stk + simulate(NW.q1f, 1000)
NW.q1mc <- fit6s.mc
NW.q1mcmc <- as.mcmc(NW.q1mc)
NW.q1smc <- NW.stk + NW.q1mc
  
NW.q2f <- fit6c
NW.q2r <- residuals(NW.q2f, NW.stk, NW.idx)
NW.q2s <- NW.stk + simulate(NW.q2f, 1000)
NW.q2mc <- fit6c.mc
NW.q2mcmc <- as.mcmc(NW.q2mc)
NW.q2smc <- NW.stk + NW.q2mc

NW.q3f <- fit6sw
NW.q3r <- residuals(NW.q3f, NW.stk, NW.idxw)
NW.q3s <- NW.stk + simulate(NW.q3f, 1000)
NW.q3mc <- fit6sw.mc
NW.q3mcmc <- as.mcmc(NW.q3mc)
NW.q3smc <- NW.stk + NW.q3mc

NW.q4f <- fit6cw
NW.q4r <- residuals(NW.q2f, NW.stk, NW.idxw)
NW.q4s <- NW.stk + simulate(NW.q4f, 1000)
NW.q4mc <- fit6cw.mc
NW.q4mcmc <- as.mcmc(NW.q4mc)
NW.q4smc <- NW.stk + NW.q4mc

NW.q5f <- fit6co
NW.q5r <- residuals(NW.q5f, NW.stk, NW.idx[1])
NW.q5s <- NW.stk + simulate(NW.q5f, 1000)
NW.q5mc <- fit6co.mc
NW.q5mcmc <- as.mcmc(NW.q5mc)
NW.q5smc <- NW.stk + NW.q5mc

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

# playing with mcmc outsputs:e.g. fit6c

qmod <- list( ~ 1, ~ 1)
fmod <- ~ s(age, k=5) + te(age, year, k=c(3,12))
m0 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=1))
m1 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250))
m2 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=500000,mcsave=250))
m3 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=500000,mcsave=1000))
# mcprobe between 0.00001 and 0.499. default=0.05. larger values fatter tails in the proposals. 
m4 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=50000,mcsave=250,mcprobe=0.1)) 
m5 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=50000,mcsave=250,mcprobe=0.2))
m6 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=50000,mcsave=250,mcprobe=0.3))
m7 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=50000,mcsave=250,mcprobe=0.4))
m8 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=50000,mcsave=250,mcprobe=0.49))
# mcrb (rescaled bounded) must be an integer between 1 and 9. default.
# it alters the covariance matrix used to propose new parameter sets in the Metropolis-Hastings algorithm
# lower values lead to a bigger reduction in correlation
m9 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=50000,mcsave=250,mcprobe=0.4, mcrb=1))
m10 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=50000,mcsave=250,mcprobe=0.4, mcrb=3))
m11 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=50000,mcsave=250,mcprobe=0.4, mcrb=5))
m12 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=50000,mcsave=250,mcprobe=0.4, mcrb=7))
m13 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=50000,mcsave=250,mcprobe=0.4, mcrb=9))


h1 <- a4aSCA(NW.stk, NW.idx, fmodel=fmod, qmodel=qmod, fit='MCMC',mcmc=SCAMCMC(mcmc=100, hybrid=T, hynstep=10, hyeps=0.1))



cbind(fitSumm(m0), fitSumm(m1), fitSumm(m2), fitSumm(m3), fitSumm(m4), 
      fitSumm(m5), fitSumm(m6), fitSumm(m7), fitSumm(m8), fitSumm(m5),
      fitSumm(m6), fitSumm(m7), fitSumm(m8), fitSumm(m9), fitSumm(m10),
      fitSumm(m11), fitSumm(m12), fitSumm(m13)   )#, fitSumm(m14), fitSumm(m5),)

plot(FLStocks(case0=NW.stk+m0, case1=NW.stk+m1, case2=NW.stk+m2, case3=NW.stk+m3))
plot(FLStocks(case2=NW.stk+m2, case4=NW.stk+m4, case5=NW.stk+m5, case6=NW.stk+m6, case7=NW.stk+m7, case8=NW.stk+m8))


pdf("m0.pdf")
out <- as.mcmc(m0)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m1.pdf")
out <- as.mcmc(m1)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m2.pdf")
out <- as.mcmc(m2)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m3.pdf")
out <- as.mcmc(m3)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()


pdf("m4.pdf")
out <- as.mcmc(m4)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m5.pdf")
out <- as.mcmc(m5)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m6.pdf")
out <- as.mcmc(m6)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m7.pdf")
out <- as.mcmc(m7)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m8.pdf")
out <- as.mcmc(m8)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m9.pdf")
out <- as.mcmc(m9)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m10.pdf")
out <- as.mcmc(m10)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m11.pdf")
out <- as.mcmc(m11)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m12.pdf")
out <- as.mcmc(m12)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

pdf("m13.pdf")
out <- as.mcmc(m13)
plot(out)
autocorr.plot(out)
levelplot(cor(out), col.regions=tim.colors(100), cuts=100, at=seq(-1, 1, length.out=100))
dev.off()

save.image("main.RData")