
   
load("C:\\Leire\\Sardina\\Ispra Whorkshop\\Ernesto\\BBdata.RData")

##original data 8+ group##############
############## plusgroup 6############
library(copula);library(triangle);library(mgcv);library(splines);library(plyr);library(ggplot2);library(knitr)

library(FLCore);library(FLEDA);library(FLash);library(ggplotFL);library(FLBRP);library(FLa4a)

library(FLXSA); library(FLEDA);library(FLa4a);library(diagram)

BB.stk<-setPlusGroup(BB.stk,6)
BB.idx[[1]]<-FLIndex(index=setPlusGroup(index(BB.idx[[1]]),6,by='sum'))
range(BB.idx[[1]])[c('min', 'max','startf', 'endf')]<-c(1,6,4.8/12,6/12) #noiz egin den: 2013 April 24st to June 3th
range(BB.idx[[2]])[c('min', 'max','startf', 'endf')]<-c(1,6,5.1/12,5.8/12)

mat(BB.stk)[2,]<-1 #bayesianarekin bat etortzeko

m.spwn(BB.stk)<-0
harvest.spwn(BB.stk)<-0
range(BB.stk)[c('minfbar','maxfbar')]<-c(2,5)

################################################################################################################
  

#Model1

#index.var(BB.idx[[1]])<-0.5

#fmodel<-~te(replace(age,age>5,5), year, k = c(6,7))

#qmodel<-list(~ s(age, k=6)+factor(as.numeric(year %in% c(2003, 2006:2007))),~1)

 
#################################
########FIT 4###################

###FINAL MODEL1##################
index.var(BB.idx[[1]])<-0.5
fmodel<-~te(replace(age,age>5,5), year, k = c(6,7))
qmodel<-list(~ s(age, k=3)+factor(as.numeric(year %in% c(2003, 2006:2007))),~1)
BB.fit4.1 <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel)

BB.fit4.1mc <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3,mcrb=7))


BB.fit4<-BB.fit4.1

##plots
plot(residuals(BB.fit4, BB.stk, BB.idx))
bubbles(residuals(BB.fit4, BB.stk, BB.idx))

stk <- BB.stk + BB.fit4
plot(stk)

plot(BB.fit4,BB.stk)
plot(BB.fit4,BB.idx[1])


a<-wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality")#, screen = list(x = -90, y=-45))
a1<-a

b<-wireframe(data ~ age + year, data = as.data.frame(predict(BB.fit4)$qmodel[[1]]), drape = TRUE, main="q acoustic")#, screen = list(x = -90, y=-45))

# arrange the two plots horizonatal
print(a, position=c(-0.5, 0, 1, 0.7), more=TRUE)
print(b, position=c(0.5, 0, 1, 0.7))

fitSumm(BB.fit4)
  

#Model2

#index.var(BB.idx[[1]])<-1

#fmodel<-~te(age, year, k = c(6,7))

#qmodel<-list(~ s(age, k=6),~1)

 

###FINAL MODEL2################## BB.q1
index.var(BB.idx[[1]])<-1
fmodel<-~te(age, year, k = c(6,7))
qmodel<-list(~ s(age, k=3),~1)
BB.fit4.2 <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel)

BB.fit4.2mc <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3,mcrb=7))

BB.q1f<-BB.fit4.2
BB.q1mc<-BB.fit4.2mc
BB.q1r<-residuals(BB.fit4.2, BB.stk, BB.idx)
BB.q1s<-BB.stk + BB.fit4.2
BB.q1smc<-BB.stk+BB.fit4.2mc
BB.q1mcmc<-as.mcmc(BB.fit4.2mc)


BB.fit4<-BB.fit4.2

##plots
plot(residuals(BB.fit4, BB.stk, BB.idx))
bubbles(residuals(BB.fit4, BB.stk, BB.idx))

stk <- BB.stk + BB.fit4
plot(stk)

plot(BB.fit4,BB.stk)
plot(BB.fit4,BB.idx[1])


a<-wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality")#, screen = list(x = -90, y=-45))
a2<-a

b<-wireframe(data ~ age + year, data = as.data.frame(predict(BB.fit4)$qmodel[[1]]), drape = TRUE, main="q acoustic")#, screen = list(x = -90, y=-45))

# arrange the two plots horizonatal
print(a, position=c(-0.5, 0, 1, 0.7), more=TRUE)
print(b, position=c(0.5, 0, 1, 0.7))
  


#Model3

#index.var(BB.idx[[1]])<-1

#fmodel<-~te(age, year, k = c(6,7))

#qmodel<-list(~ 1,~1)

 
###FINAL MODEL3################## BB.q2
index.var(BB.idx[[1]])<-1
fmodel<-~te(age, year, k = c(6,7))
qmodel<-list(~ 1,~1)
BB.fit4.3 <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel)

BB.fit4.3mc <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3,mcrb=7))


BB.q2f<-BB.fit4.3
BB.q2mc<-BB.fit4.3mc
BB.q2r<-residuals(BB.fit4.3, BB.stk, BB.idx)
BB.q2s<-BB.stk + BB.fit4.3
BB.q2smc<-BB.stk+BB.fit4.3mc
BB.q2mcmc<-as.mcmc(BB.fit4.3mc)



BB.fit4<-BB.fit4.3

##plots
plot(residuals(BB.fit4, BB.stk, BB.idx))
bubbles(residuals(BB.fit4, BB.stk, BB.idx))

stk <- BB.stk + BB.fit4
plot(stk)

plot(BB.fit4,BB.stk)
plot(BB.fit4,BB.idx[1])

a<-wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality")#, screen = list(x = -90, y=-45))
a3<-a

b<-wireframe(data ~ age + year, data = as.data.frame(predict(BB.fit4)$qmodel[[1]]), drape = TRUE, main="q acoustic")#, screen = list(x = -90, y=-45))

# arrange the two plots horizonatal
print(a, position=c(-0.5, 0, 1, 0.7), more=TRUE)
print(b, position=c(0.5, 0, 1, 0.7))
  

#Model4

#index.var(BB.idx[[1]])<-0.5

#fmodel<-~te(age, year, k = c(6,7))

#qmodel<-list(~ s(age, k=6),~1)

 

###FINAL MODEL4################## BB.q3
index.var(BB.idx[[1]])<-0.5
fmodel<-~te(age, year, k = c(6,7))
qmodel<-list(~ s(age, k=3),~1)
BB.fit4.4 <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel)

BB.fit4.4mc <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3,mcrb=7))

BB.q3f<-BB.fit4.4
BB.q3mc<-BB.fit4.4mc
BB.q3r<-residuals(BB.fit4.4, BB.stk, BB.idx)
BB.q3s<-BB.stk + BB.fit4.4
BB.q3smc<-BB.stk+BB.fit4.4mc
BB.q3mcmc<-as.mcmc(BB.fit4.4mc)

BB.fit4<-BB.fit4.4



##plots
plot(residuals(BB.fit4, BB.stk, BB.idx))
bubbles(residuals(BB.fit4, BB.stk, BB.idx))

stk <- BB.stk + BB.fit4
plot(stk)

plot(BB.fit4,BB.stk)
plot(BB.fit4,BB.idx[1])


a<-wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality")#, screen = list(x = -90, y=-45))
a2<-a

b<-wireframe(data ~ age + year, data = as.data.frame(predict(BB.fit4)$qmodel[[1]]), drape = TRUE, main="q acoustic")#, screen = list(x = -90, y=-45))

# arrange the two plots horizonatal
print(a, position=c(-0.5, 0, 1, 0.7), more=TRUE)
print(b, position=c(0.5, 0, 1, 0.7))
  


#Model5

#index.var(BB.idx[[1]])<-0.5

#fmodel<-~te(age, year, k = c(6,7))

#qmodel<-list(~ 1,~1)

 
###FINAL MODEL5################## BB.q4
index.var(BB.idx[[1]])<-0.5
fmodel<-~te(age, year, k = c(6,7))
qmodel<-list(~ 1,~1)
BB.fit4.5 <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel)

BB.fit4.5mc <- a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3,mcrb=7))


BB.q4f<-BB.fit4.5
BB.q4mc<-BB.fit4.5mc
BB.q4r<-residuals(BB.fit4.5, BB.stk, BB.idx)
BB.q4s<-BB.stk + BB.fit4.5
BB.q4smc<-BB.stk+BB.fit4.5mc
BB.q4mcmc<-as.mcmc(BB.fit4.5mc)

BB.fit4<-BB.fit4.5

##plots
plot(residuals(BB.fit4, BB.stk, BB.idx))
bubbles(residuals(BB.fit4, BB.stk, BB.idx))

stk <- BB.stk + BB.fit4
plot(stk)

plot(BB.fit4,BB.stk)
plot(BB.fit4,BB.idx[1])

a<-wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality")#, screen = list(x = -90, y=-45))
a3<-a

b<-wireframe(data ~ age + year, data = as.data.frame(predict(BB.fit4)$qmodel[[1]]), drape = TRUE, main="q acoustic")#, screen = list(x = -90, y=-45))

# arrange the two plots horizonatal
print(a, position=c(-0.5, 0, 1, 0.7), more=TRUE)
print(b, position=c(0.5, 0, 1, 0.7))
  


#ALL


 
####all

plot(FLStocks(selected=BB.stk+BB.fit4.1,f_te.q_s=BB.stk+BB.fit4.2,f_te.q_1=BB.stk+BB.fit4.3,W_f_te.q_s=BB.stk+BB.fit4.4,W_f_te.q_1=BB.stk+BB.fit4.5,m6=BB.stk+BB.fit4.6,m7=BB.stk+BB.fit4.7))

plot(FLStocks(selected=BB.stk+BB.fit4.1mc,f_te.q_s=BB.stk+BB.fit4.2mc,f_te.q_1=BB.stk+BB.fit4.3mc,W_f_te.q_s=BB.stk+BB.fit4.4mc,W_f_te.q_1=BB.stk+BB.fit4.5mc,m6=BB.stk+BB.fit4.6mc,m7=BB.stk+BB.fit4.7mc))

  


#MCMC Correlation  for Model 1

mcmcFit<-as.mcmc(BB.fit4.1mc)
levelplot(cor(BB.q4mcmc),col.regions=rainbow(1000)[600:1],cuts=100,at=seq(-1,1, length.out=600))
  

