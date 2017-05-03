################################################
# Stock assessment of the full area (BB+NW+S)
# - Data compilation 
# - Model runs
################################################


library(FLa4a)

load("../BB/BBdata.RData")
load("../IB/IB.Rdata")
load("../NW/NW.Rdata")
load("../S/S.RData")

#### Set plusgroup for BoB sardine ####
BB.stkpg <- setPlusGroup(BB.stk, 6)


#### Trimming #########################
yearsel=2002:2014
BB.temp<-trim(BB.stkpg, year=yearsel)
IB.temp<-trim(IB.stk, year=yearsel)
NW.temp<-trim(NW.stk, year=yearsel)
S.temp<-trim(S.stk, year=yearsel)


#### Building A.stk object ############
A.stk<-FLStock()
name(A.stk)<-"ALL REGIONS ATLANTO-IBERIAN SARDINE"
desc(A.stk)<-"prepared for a4a workshop - 15/12/2015"
A.stk@range<-IB.temp@range
A.stk@m.spwn<-IB.temp@m.spwn
A.stk@harvest.spwn<-IB.temp@harvest.spwn
A.stk@harvest<-BB.temp@harvest
A.stk@mat<- S.temp@mat
A.stk@m<-S.temp@m
A.stk@discards<-BB.temp@discards
A.stk@discards.wt<-BB.temp@discards.wt
A.stk@discards.n<-BB.temp@discards.n
A.stk@landings<- BB.temp@landings+S.temp@landings*1000+NW.temp@landings*1000
units(A.stk@landings)<-"tons"
A.stk@landings.n<- BB.temp@landings.n+S.temp@landings.n*1000+NW.temp@landings.n*1000
units(A.stk@landings.n)<-"thousands"
A.stk@catch<-A.stk@landings
A.stk@catch.n<- A.stk@landings.n
A.stk@landings.wt<-  (BB.temp@catch.wt*BB.temp@landings.n + NW.temp@catch.wt*NW.temp@landings.n*1000+ S.temp@catch.wt*S.temp@landings.n*1000)/A.stk@landings.n
units(A.stk@landings.wt)<-"kg"
A.stk@catch.wt<-A.stk@landings.wt
A.stk@stock<- S.temp@stock
A.stk@stock.n<- S.temp@stock.n
A.stk@stock.wt<-  (BB.temp@stock.wt*BB.temp@landings.n + NW.temp@stock.wt*NW.temp@landings.n*1000+ S.temp@stock.wt*S.temp@landings.n*1000)/A.stk@landings.n
units(A.stk@stock.wt)<-"kg"


#### Building A.idx object ############
BB1<-trim(BB.idx[[1]],year=yearsel,age=1:6)
BB2<-trim(BB.idx[[2]],year=yearsel,age=1:6)
IB1<-trim(IB.idx[[1]],year=yearsel,age=1:6)
IB2<-trim(IB.idx[[2]],year=yearsel,age=1:6)
NW1<-trim(NW.idx[[1]],year=yearsel,age=1:6)
NW2<-trim(NW.idx[[2]],year=yearsel,age=1:6)
S1 <-trim( S.idx[[1]],year=yearsel,age=1:6)
S2 <-trim( S.idx[[2]],year=yearsel,age=1:6)
desc(BB2)<-"BB Egg"
desc(BB1)<-"BB Acoustic"
desc(IB1)<-"IB Acoustic"
desc(IB2)<-"IB DEPM"
desc(NW1)<-"NW Acoustic"
desc(NW2)<-"NW DEPM"
desc(S1)<-"S Acoustic"
desc(S2)<-"S DEPM"
name(BB2)<-"BB Egg"
name(BB1)<-"BB Acoustic"
name(IB1)<-"IB Acoustic"
name(IB2)<-"IB DEPM"
name(NW1)<-"NW ACoustic"
name(NW2)<-"NW DEPM"
name(S1)<-"S Acoustic"
name(S2)<-"S DEPM"


#### Seperate survey file  ############

### Separate object file ###
NW1@index[,3]<- (NW1@index[,2]+NW1@index[,4])/2
NW1@index[,11]<- (NW1@index[,10]+NW1@index[,12])/2
S1@index[,3]<- (S1@index[,2]+S1@index[,4])/2
S1@index[,11]<- (S1@index[,10]+S1@index[,12])/2

#### AC: Aggregated acoustic index ####
AC<-NW1
name(AC)<-"Aggr Acoustic"
desc(AC)<-"Aggr Acoustic"
AC@index<-BB1@index+ NW1@index*1000 + S1@index*1000
AC@range<-BB1@range
units(AC@index)<-"Thousands"
A.idx<-FLIndices(AC,BB2,NW2,S2) 

#### DP: Aggregated DEPM index ########
DP<-NW2
name(DP)<-"Aggr DEPM"
desc(DP)<-"Aggr DEPM"
DP@index<-NW2@index*1000 + S2@index*1000
DP@range<-NW2@range
units(DP@index)<-"Thousands"
DP@index[,2]<-(DP@index[,1]*2+DP@index[,4]*1)/3
#DP@index[,3]<-(DP@index[,1]*1+DP@index[,4]*2)/3
#DP@index[,5]<-(DP@index[,4]*2+DP@index[,7]*1)/3
#DP@index[,6]<-(DP@index[,4]*1+DP@index[,7]*2)/3
#DP@index[,8]<-(DP@index[,7]*2+DP@index[,10]*1)/3
#DP@index[,9]<-(DP@index[,7]*1+DP@index[,10]*2)/3
#DP@index[,11]<-(DP@index[,10]*2+DP@index[,13]*1)/3
#DP@index[,12]<-(DP@index[,10]*1+DP@index[,13]*2)/3


A.idx<-FLIndices(AC,BB2,DP) 		### Aggr Acoustic + DEPM
A.idx2s<-FLIndices(BB1, NW1, S1,DP) 	### Aggr DEPM but separate acoustic
A.idx5<-FLIndices(AC) 			### Acoustic only






#### Model Runs ##########################
temp.idx<-A.idx
index.var(temp.idx[[1]]) ### Make sure we have no weighting in the survey: look for NA ###

fmodel  <- ~ te(age, year, k=c(3,7)) + s(age, k=5)

q1 	<- list(~s(age,k=4),~1,~1)				# Standard run

q2 	<- list(~factor(replace(age,age>1,2)), ~year,~1)	# Seperate catchability between age 1 and older ages
q2extra <- list(~1, ~year,~1)					# linear year effect on BB egg survey
q2split <- list(~1, ~1, ~1, ~1)					# Constant catchability accross all (disaggregated acoustic surveys and sum of DEPM)
q2splity<- list(~year, ~year, ~1, ~1)				# linear year effect on BB and NW acoustic survey
q5 	<- list(~1)						# just a single aggregated acoustic survey

A.q1f 		<- a4aSCA(A.stk, A.idx,fmodel=fmodel,qmodel=q1)
A.q2f 		<- a4aSCA(A.stk, A.idx,fmodel=fmodel,qmodel=q2)
A.q2extraf 	<- a4aSCA(A.stk, A.idx,fmodel=fmodel,qmodel=q2extra)
A.q2splitf 	<- a4aSCA(A.stk, A.idx2s,fmodel=fmodel,qmodel=q2split)
A.q2splityf 	<- a4aSCA(A.stk, A.idx2s,fmodel=fmodel,qmodel=q2splity)
A.q5f 		<- a4aSCA(A.stk, A.idx5,fmodel=fmodel,qmodel=q5)

A.q1mc 		<- a4aSCA(A.stk, A.idx,fmodel=fmodel,qmodel=q1,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3))
A.q2mc 		<- a4aSCA(A.stk, A.idx,fmodel=fmodel,qmodel=q2,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3))
A.q2extramc 	<- a4aSCA(A.stk, A.idx,fmodel=fmodel,qmodel=q2extra,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3))
A.q2splitmc 	<- a4aSCA(A.stk, A.idx2s,fmodel=fmodel,qmodel=q2split,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3))
A.q2splitymc 	<- a4aSCA(A.stk, A.idx2s,fmodel=fmodel,qmodel=q2splity,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3))
A.q5mc 		<- a4aSCA(A.stk, A.idx5,fmodel=fmodel,qmodel=q5,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3))

A.q1r 		<- residuals(A.q1f, A.stk, A.idx)
A.q2r 		<- residuals(A.q2f, A.stk, A.idx)
A.q2extrar 	<- residuals(A.q2extraf, A.stk, A.idx)
A.q2splitr 	<- residuals(A.q2splitf, A.stk, A.idx2s)
A.q2splityr 	<- residuals(A.q2splityf, A.stk, A.idx2s)
A.q5r 		<- residuals(A.q5f, A.stk, A.idx5)

A.q1mcmc	<- as.mcmc(A.q1mc)
A.q2mcmc	<- as.mcmc(A.q2mc)
A.q2extramcmc	<- as.mcmc(A.q2extramc)
A.q2splitmcmc	<- as.mcmc(A.q2splitmc)
A.q2splitymcmc	<- as.mcmc(A.q2splitymc)
A.q5mcmc	<- as.mcmc(A.q5mc)


stk <- A.stk + A.q1f
A.q1s <- stk + simulate(A.q1f, 500)
stk <- A.stk + A.q2f
A.q2s <- stk + simulate(A.q2f, 500)
stk <- A.stk + A.q2extraf
A.q2extras <- stk + simulate(A.q2extraf, 500)
stk <- A.stk + A.q2splitf
A.q2splits <- stk + simulate(A.q2splitf, 500)
stk <- A.stk + A.q2splityf
A.q2splitys <- stk + simulate(A.q2splityf, 500)
A.q2splitysmc <- A.stk + A.q2splitymc
stk <- A.stk + A.q5f
A.q5s <- stk + simulate(A.q5f, 500)
A.q5smc <- A.stk + A.q5mc


A.q1smc <- A.stk + A.q1mc 
A.q2smc <- A.stk + A.q2mc
A.q2extrasmc <- A.stk + A.q2extramc
A.q2splitsmc <- A.stk + A.q2splitmc




############
index.var(A.idx[[1]])<-0.5  ##### Change the weight of the aggregated acoustic survey ###
############


A.q3f 		<- a4aSCA(A.stk, A.idx,fmodel=fmodel,qmodel=q1) # same q than q1 run but survey weight changed
A.q4f 		<- a4aSCA(A.stk, A.idx,fmodel=fmodel,qmodel=q2) # same q than q2 run but survey weight changed 

A.q3mc 		<- a4aSCA(A.stk, A.idx,fmodel=fmodel,qmodel=q1,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3))
A.q4mc 		<- a4aSCA(A.stk, A.idx,fmodel=fmodel,qmodel=q2,fit='MCMC',mcmc=SCAMCMC(mcmc=12500,mcsave=250,mcprobe=0.3))

A.q3r 		<- residuals(A.q3f, A.stk, A.idx)
A.q4r 		<- residuals(A.q4f, A.stk, A.idx)

A.q3mcmc	<- as.mcmc(A.q3mc)
A.q4mcmc	<- as.mcmc(A.q4mc)

stk <- A.stk + A.q3f
A.q3s <- stk + simulate(A.q3f, 500)
stk <- A.stk + A.q4f
A.q4s <- stk + simulate(A.q4f, 500)

A.q3smc <- A.stk + A.q3mc 
A.q4smc <- A.stk + A.q4mc


save(file="Amodels.Rdata",A.idx, A.idx2s, A.idx5, A.stk, 
		A.q1f,A.q2f,A.q3f,A.q4f,A.q2extraf,A.q2splitf,A.q2splityf,A.q5f,
		A.q1r,A.q2r,A.q3r,A.q4r,A.q2extrar,A.q2splitr,A.q2splityr,A.q5r,
		A.q1s,A.q2s,A.q3s,A.q4s,A.q2extras,A.q2splits,A.q2splitys,A.q5s,
		A.q1mc,A.q2mc,A.q3mc,A.q4mc,A.q2extramc,A.q2splitmc,A.q2splitymc,A.q5mc,
		A.q1mcmc,A.q2mcmc,A.q3mcmc,A.q4mcmc,A.q2extramcmc,A.q2splitmcmc,A.q2splitymcmc,A.q5mcmc,
		A.q1smc, A.q2smc, A.q3smc, A.q4smc,A.q2extrasmc,A.q2splitsmc, A.q2splitysmc, A.q5smc)



