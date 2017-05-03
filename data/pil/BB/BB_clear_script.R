load("C:\\Leire\\Sardina\\Ispra Whorkshop\\Ernesto\\BBdata.RData")

##original data 8+ group##############
############## plusgroup 6############

BB.stk<-setPlusGroup(BB.stk,6)
BB.idx[[1]]<-FLIndex(index=setPlusGroup(index(BB.idx[[1]]),6,by='sum'))
range(BB.idx[[1]])[c('min', 'max','startf', 'endf')]<-c(1,6,4.8/12,6/12) #noiz egin den: 2013 April 24st to June 3th
range(BB.idx[[2]])[c('min', 'max','startf', 'endf')]<-c(1,6,5.1/12,5.8/12)

mat(BB.stk)[2,]<-1 #bayesianarekin bat etortzeko

m.spwn(BB.stk)<-0
harvest.spwn(BB.stk)<-0
range(BB.stk)[c('minfbar','maxfbar')]<-c(2,5)


###FINAL MODEL2################## BB.q1
index.var(BB.idx[[1]])<-1
fmodel<-~te(age, year, k = c(6,7))
qmodel<-list(~ s(age, k=3),~1)

BB.q1f<-a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel)
BB.q1mc<-a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3,mcrb=7))
BB.q1r<-residuals(BB.q1f, BB.stk, BB.idx)
BB.q1s<-BB.stk + BB.q1f
BB.q1smc<-BB.stk+BB.q1mc
BB.q1mcmc<-as.mcmc(BB.q1mc)


BB.fit4<-BB.q1f
BB.fit4mc<-BB.q1mc


##plots
plot(residuals(BB.fit4, BB.stk, BB.idx))
bubbles(residuals(BB.fit4, BB.stk, BB.idx))

stk <- BB.stk + BB.fit4mc
plot(stk)

plot(BB.fit4,BB.stk)
plot(BB.fit4,BB.idx[1])


a<-wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality")#, screen = list(x = -90, y=-45))


b<-wireframe(data ~ age + year, data = as.data.frame(predict(BB.fit4)$qmodel[[1]]), drape = TRUE, main="q acoustic")#, screen = list(x = -90, y=-45))

# arrange the two plots horizontal
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

BB.q2f<-a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel)
BB.q2mc<-a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3,mcrb=7))
BB.q2r<-residuals(BB.q2f, BB.stk, BB.idx)
BB.q2s<-BB.stk + BB.q2f
BB.q2smc<-BB.stk+BB.q2mc
BB.q2mcmc<-as.mcmc(BB.q2mc)


BB.fit4<-BB.q2f
BB.fit4mc<-BB.q2mc

##plots
plot(residuals(BB.fit4, BB.stk, BB.idx))
bubbles(residuals(BB.fit4, BB.stk, BB.idx))

stk <- BB.stk + BB.fit4mc
plot(stk)

plot(BB.fit4,BB.stk)
plot(BB.fit4,BB.idx[1])

a<-wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality")#, screen = list(x = -90, y=-45))


b<-wireframe(data ~ age + year, data = as.data.frame(predict(BB.fit4)$qmodel[[1]]), drape = TRUE, main="q acoustic")#, screen = list(x = -90, y=-45))

# arrange the two plots horizontal
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

BB.q3f<-a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel)
BB.q3mc<-a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3,mcrb=7))
BB.q3r<-residuals(BB.q3f, BB.stk, BB.idx)
BB.q3s<-BB.stk + BB.q3f
BB.q3smc<-BB.stk+BB.q3mc
BB.q3mcmc<-as.mcmc(BB.q3mc)


BB.fit4<-BB.q3f
BB.fit4mc<-BB.q3mc



##plots
plot(residuals(BB.fit4, BB.stk, BB.idx))
bubbles(residuals(BB.fit4, BB.stk, BB.idx))

stk <- BB.stk + BB.fit4
plot(stk)

plot(BB.fit4,BB.stk)
plot(BB.fit4,BB.idx[1])


a<-wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality")#, screen = list(x = -90, y=-45))


b<-wireframe(data ~ age + year, data = as.data.frame(predict(BB.fit4)$qmodel[[1]]), drape = TRUE, main="q acoustic")#, screen = list(x = -90, y=-45))

# arrange the two plots horizontal
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

BB.q4f<-a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel)
BB.q4mc<-a4aSCA(BB.stk, BB.idx,fmodel=fmodel,qmodel=qmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3,mcrb=7))
BB.q4r<-residuals(BB.q4f, BB.stk, BB.idx)
BB.q4s<-BB.stk + BB.q4f
BB.q4smc<-BB.stk+BB.q4mc
BB.q4mcmc<-as.mcmc(BB.q4mc)


BB.fit4<-BB.q4f
BB.fit4mc<-BB.q4mc

##plots
plot(residuals(BB.fit4, BB.stk, BB.idx))
bubbles(residuals(BB.fit4, BB.stk, BB.idx))

stk <- BB.stk + BB.fit4mc
plot(stk)

plot(BB.fit4,BB.stk)
plot(BB.fit4,BB.idx[1])

a<-wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality")#, screen = list(x = -90, y=-45))
b<-wireframe(data ~ age + year, data = as.data.frame(predict(BB.fit4)$qmodel[[1]]), drape = TRUE, main="q acoustic")#, screen = list(x = -90, y=-45))

# arrange the two plots horizontal
print(a, position=c(-0.5, 0, 1, 0.7), more=TRUE)
print(b, position=c(0.5, 0, 1, 0.7))



###FINAL MODEL6################## BB.q5 (q2 without depm)
index.var(BB.idx[[1]])<-1
fmodel<-~te(age, year, k = c(6,7))
qmodel<-list(~ 1)

BB.q5f<-a4aSCA(BB.stk, BB.idx[1],fmodel=fmodel,qmodel=qmodel)
BB.q5mc<-a4aSCA(BB.stk, BB.idx[1],fmodel=fmodel,qmodel=qmodel,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3,mcrb=7))
BB.q5r<-residuals(BB.q5f, BB.stk, BB.idx[1])
BB.q5s<-BB.stk + BB.q5f
BB.q5smc<-BB.stk+BB.q5mc
BB.q5mcmc<-as.mcmc(BB.q5mc)

BB.fit4<-BB.q5f
BB.fit4mc<-BB.q5mc

##plots
plot(residuals(BB.fit4, BB.stk, BB.idx[1]))
bubbles(residuals(BB.fit4, BB.stk, BB.idx[1]))

stk <- BB.stk + BB.fit4mc
plot(stk)

plot(BB.fit4,BB.stk)
plot(BB.fit4,BB.idx[1])

a<-wireframe(data ~ age + year, data = as.data.frame(harvest(stk)), drape = TRUE, main="Fishing mortality")#, screen = list(x = -90, y=-45))
b<-wireframe(data ~ age + year, data = as.data.frame(predict(BB.fit4)$qmodel[[1]]), drape = TRUE, main="q acoustic")#, screen = list(x = -90, y=-45))

# arrange the two plots horizontal
print(a, position=c(-0.5, 0, 1, 0.7), more=TRUE)
print(b, position=c(0.5, 0, 1, 0.7))



save(list=c("BB.stk","BB.idx", paste("BB.", rep(c("q1","q2","q3","q4","q5"), each=6), rep(c("f","r","s","mc","mcmc","smc"),4), sep="")),file='C:\\Leire\\Sardina\\Ispra Whorkshop\\clear_object_results.RData')

