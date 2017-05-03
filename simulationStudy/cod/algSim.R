library(FLBRP)
library(FLa4a)
library(ggplotFL)
source("funs.R")
library(parallel)

#====================================================================
# Simulation study spatial heterogeneity
#====================================================================

#--------------------------------------------------------------------
# Simulate # pars loosely based on fishbase and ICES
#--------------------------------------------------------------------

linf <- 140
k <- 0.125
t0 <- 0
a1 <- 2
a50 <- 3.5
stp <- seq(0.7,0.95, 0.05)
sl <- 2
sr <- 100
v <- 1000
cvf <- 0.4
cvr <- 0.4
cvq <- 0.4
cvc <- 0.4
ny <- 30
mc <- 21
maxage <- 15
it <- 250
seed <- 1133
set.seed(seed)
srr <- FLQuant(rlnorm(ny*it, 0, log(cvr^2 + 1)), dimnames=list(year=1:ny, iter=1:it))

oms <- mclapply(split(stp, stp), function(x, srModel="bevholt"){
	brp <- lh(gislasim(linf=linf, k=k, t0=t0, s=x, a1=a1, a50=a50, sl=sl, sr=sr, v=v), range=c(min = 1, max = maxage, minfbar = 1, maxfbar = 5, plusgroup = maxage), sr=srModel)
	#m(brp)[] <- c(1.25, 0.9, 0.3, rep(0.2, dims(brp)$age-3))
	Fmsy <- refpts(brp)["msy","harvest"]
	Fc <- Fmsy*1.5
	Ftrg <- c(seq(0, Fc, len=ny*2/5), rep(Fc, ny*2/5), seq(Fc, Fmsy, len=ny/5))
	set.seed(seed)
	Ftrg <- Ftrg*rlnorm(length(Ftrg), 0, log(cvf^2+1))
	trg <- fwdControl(data.frame(year=c(1:ny+1), quantity=rep('f', ny), val=Ftrg))
	ex.stk <- as(brp, "FLStock")[,1:(ny+1)]
	catch.n(ex.stk) <- catch.n(ex.stk)[,,,,,rep(1,it)]
	# S/R
	ex.sr <- FLSR(model=srModel)
 	params(ex.sr) <- FLPar(abPars(srModel, c(refpts(brp)["virgin","ssb"]/refpts(brp)["virgin","rec"]), v=v, s=x))
	stk <- fwd(ex.stk, ctrl=trg, sr=ex.sr, sr.residuals=srr)
	list(brp=brp, stk=window(stk, ny/5, ny))
}, mc.cores=mc)

#--------------------------------------------------------------------
# no diffusion
#--------------------------------------------------------------------

# pop split
age <- as.numeric(dimnames(oms[[1]]$stk@catch.n)$age)
vps <- gle(age, b=1, q=6, m=0, nu=0.1)
vps[] <- 0.5

# sub-units
oms.su <- expand.grid(stk1=names(oms), stk2=names(oms))
nodup <- oms.su
nodup$id <- as.numeric(ac(nodup$stk1))*as.numeric(ac(nodup$stk2))
nodup <- nodup[order(nodup$id),]
nodup <- nodup[nodup$id[-nrow(nodup)]-nodup$id[-1]!=0,]
nodup$id <- paste(nodup$stk1, nodup$stk2, sep=".")

oms.su <- mclapply(split(oms.su, oms.su), function(x){
	stka <- oms[[ac(x$stk1)]]$stk
	stkb <- oms[[ac(x$stk2)]]$stk
	set.seed(seed)
	cres <- rlnorm(prod(dim(catch.n(stka))), 0, log(cvc^2 + 1))
	catch.n(stka) <- catch.n(stka)*cres 
	catch.n(stkb) <- catch.n(stkb)*cres 
	stk <- merge(stka, stkb)
	popspl <- catch.n(stk)
	popspl[] <- vps
	stkc1 <- stk
	stock.n(stkc1) <- popspl*(stock.n(stk))
	catch.n(stkc1) <- catch.n(stk)/stock.n(stk)*stock.n(stkc1)
	catch(stkc1) <- computeCatch(stkc1)
	stkc2 <- stk
	stock.n(stkc2) <- (1-popspl)*(stock.n(stk))
	catch.n(stkc2) <- catch.n(stk)/stock.n(stk)*stock.n(stkc2)
	catch(stkc2) <- computeCatch(stkc2)
	set.seed(seed)
	idxa <- getIdx(stka, cv=cvq)
	idxa <- setPG(idxa, getPG(stka))
	set.seed(seed)
	idxb <- getIdx(stkb, cv=cvq)
	idxb <- setPG(idxb, getPG(stkb))
	set.seed(seed)
	idx <- getIdx(stk, cv=cvq)
	idx <- setPG(idx, getPG(stk))
	set.seed(seed)
	idxc1 <- getIdx(stkc1, cv=cvq)
	idxc1 <- setPG(idxc1, getPG(stkc1))
	set.seed(seed)
	idxc2 <- getIdx(stkc2, cv=cvq)
	idxc2 <- setPG(idxc2, max(getPG(stkc2), 5))
	list(stk=stk, idx=idx, stka=stka, idxa=idxa, stkb=stkb, idxb=idxb, stkc1=stkc1, idxc1=idxc1, stkc2=stkc2, idxc2=idxc2)
}, mc.cores=mc)

oms.su <- oms.su[nodup$id]

#--------------------------------------------------------------------
# fit
#--------------------------------------------------------------------

fits <- mclapply(oms.su, function(x){
cat("o ") 
	fmod <- ~te(age, year, k = c(10, 10), bs = "tp") + s(year, by=as.numeric(age==1), k=5)
	qmod <- list(~s(age, k=11))
	x$fit <- a4aSCA(x[["stk"]], FLIndices(a=x[["idx"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x$fita <- a4aSCA(x[["stka"]], FLIndices(a=x[["idxa"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x$fitb <- a4aSCA(x[["stkb"]], FLIndices(a=x[["idxb"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	qmod <- list(~factor(age))
	x$fitc1 <- a4aSCA(x[["stkc1"]], FLIndices(a=x[["idxc1"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x$fitc2 <- a4aSCA(x[["stkc2"]], FLIndices(a=x[["idxc2"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	pdf(paste(format(Sys.time(), "%H:%M:%S"), "res0.pdf", sep="."))
	print(plot(residuals(x$fit, x[["stk"]], FLIndices(a=x[["idx"]]))))
	print(plot(residuals(x$fita, x[["stka"]], FLIndices(a=x[["idxa"]]))))
	print(plot(residuals(x$fitb, x[["stkb"]], FLIndices(a=x[["idxb"]]))))
	print(plot(residuals(x$fitc1, x[["stkc1"]], FLIndices(a=x[["idxc1"]]))))
	print(plot(residuals(x$fitc2, x[["stkc2"]], FLIndices(a=x[["idxc2"]]))))
	dev.off()
	x
}, mc.cores=mc)

# check failures
lst <- lapply(fits, function(x) c(anyNA(harvest(x$fit)), anyNA(harvest(x$fita)), anyNA(harvest(x$fitb)), anyNA(harvest(x$fitc1)), anyNA(harvest(x$fitc2))))
fails <- do.call("rbind", lst)

#--------------------------------------------------------------------
# Post-process
#--------------------------------------------------------------------

spatialHeterogeneityRes <- mclapply(fits, function(x){
	stks  <- with(x, FLStocks(stk=stk+fit, stka=stka+fita, stkb=stkb+fitb, stkc1=stkc1+fitc1, stkc2=stkc2+fitc2))
	stks$stksu <- merge(stks$stka, stks$stkb)
	stks$stksmx <- merge(stks$stkc1, stks$stkc2)
	gc()
	list(stks=stks)
}, mc.cores=mc)

diffRes <- mclapply(spatialHeterogeneityRes, function(x){
	ssbmseu <- with(x, median((ssb(stks$stk)-ssb(stks$stksu))^2))
	ssbmsemx <- with(x, median((ssb(stks$stk)-ssb(stks$stksmx))^2))
	recmseu <- with(x, median((rec(stks$stk)-rec(stks$stksu))^2))
	recmsemx <- with(x, median((rec(stks$stk)-rec(stks$stksmx))^2))
	data.frame(ssbmseu=ssbmseu, ssbmsemx=ssbmsemx, recmseu=recmseu, recmsemx=recmsemx)
}, mc.cores=mc)

diffRes <- cbind(do.call("rbind", diffRes), stk1=as.numeric(ac(nodup$stk1)), stk2=as.numeric(ac(nodup$stk2)))

# to be used later, can be dropped at the end
spatialHeterogeneityRes0 <- spatialHeterogeneityRes

# save objects
save(spatialHeterogeneityRes, diffRes, vps, oms.su, fits, fails, file="RData.sim0")
rm(spatialHeterogeneityRes, diffRes, oms.su, fits, fails)

#--------------------------------------------------------------------
# diffusion 0.1-0.9
#--------------------------------------------------------------------

# pop split
age <- as.numeric(dimnames(oms[[1]]$stk@catch.n)$age)
vps <- gle(age, b=1, q=6, m=0, nu=0.1)
vps[vps<0.1] <- 0.1
vps[vps>0.9] <- 0.9

# sub-units
oms.su <- expand.grid(stk1=names(oms), stk2=names(oms))
nodup <- oms.su
nodup$id <- as.numeric(ac(nodup$stk1))*as.numeric(ac(nodup$stk2))
nodup <- nodup[order(nodup$id),]
nodup <- nodup[nodup$id[-nrow(nodup)]-nodup$id[-1]!=0,]
nodup$id <- paste(nodup$stk1, nodup$stk2, sep=".")

oms.su <- mclapply(split(oms.su, oms.su), function(x){
	stka <- oms[[ac(x$stk1)]]$stk
	stkb <- oms[[ac(x$stk2)]]$stk
	set.seed(seed)
	cres <- rlnorm(prod(dim(catch.n(stka))), 0, log(cvc^2 + 1))
	catch.n(stka) <- catch.n(stka)*cres 
	catch.n(stkb) <- catch.n(stkb)*cres 
	stk <- merge(stka, stkb)
	popspl <- catch.n(stk)
	popspl[] <- vps
	stkc1 <- stk
	stock.n(stkc1) <- popspl*(stock.n(stk))
	catch.n(stkc1) <- catch.n(stk)/stock.n(stk)*stock.n(stkc1)
	catch(stkc1) <- computeCatch(stkc1)
	stkc2 <- stk
	stock.n(stkc2) <- (1-popspl)*(stock.n(stk))
	catch.n(stkc2) <- catch.n(stk)/stock.n(stk)*stock.n(stkc2)
	catch(stkc2) <- computeCatch(stkc2)
	set.seed(seed)
	idxc1 <- getIdx(stkc1, cv=cvq)
	idxc1 <- setPG(idxc1, getPG(stkc1))
	set.seed(seed)
	idxc2 <- getIdx(stkc2, cv=cvq)
	idxc2 <- setPG(idxc2, max(getPG(stkc2), 5))
	list(stkc1=stkc1, idxc1=idxc1, stkc2=stkc2, idxc2=idxc2)
}, mc.cores=mc)

oms.su <- oms.su[nodup$id]

#--------------------------------------------------------------------
# fit
#--------------------------------------------------------------------

fits <- mclapply(oms.su, function(x){
cat("o ") 
	fmod <- ~te(age, year, k = c(10, 10), bs = "tp") + s(year, by=as.numeric(age==1), k=5)
	qmod <- list(~s(age, k=11))
	x$fitc1 <- a4aSCA(x[["stkc1"]], FLIndices(a=x[["idxc1"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x$fitc2 <- a4aSCA(x[["stkc2"]], FLIndices(a=x[["idxc2"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x
}, mc.cores=mc)

# check failures
lst <- lapply(fits, function(x) c(anyNA(harvest(x$fitc1)), anyNA(harvest(x$fitc2))))
fails <- do.call("rbind", lst)

#--------------------------------------------------------------------
# Post-process
#--------------------------------------------------------------------

spatialHeterogeneityRes <- lapply(fits, function(x){
	stks  <- with(x, FLStocks(stkc1=stkc1+fitc1, stkc2=stkc2+fitc2))
	stks$stksmx <- merge(stks$stkc1, stks$stkc2)
	stks$stkp <- x$stkp
	gc()
	list(stks=stks)
})

# recover overall stocks
for(i in 1:length(spatialHeterogeneityRes)) spatialHeterogeneityRes[[i]]$stks$stk <- spatialHeterogeneityRes0[[i]]$stks$stk

diffRes <- mclapply(spatialHeterogeneityRes, function(x){
	ssbmsemx <- with(x, median((ssb(stks$stk)-ssb(stks$stksmx))^2))
	recmsemx <- with(x, median((rec(stks$stk)-rec(stks$stksmx))^2))
	data.frame(ssbmsemx=ssbmsemx, recmsemx=recmsemx)
}, mc.cores=mc)

diffRes <- cbind(do.call("rbind", diffRes), stk1=as.numeric(ac(nodup$stk1)), stk2=as.numeric(ac(nodup$stk2)))

save(spatialHeterogeneityRes, diffRes, vps, oms.su, fits, fails, file="RData.sim1")
rm(spatialHeterogeneityRes, diffRes, oms.su, fits, fails)

#--------------------------------------------------------------------
# diffusion 0.2-0.8
#--------------------------------------------------------------------

# pop split
age <- as.numeric(dimnames(oms[[1]]$stk@catch.n)$age)
vps <- gle(age, b=1, q=6, m=0, nu=0.1)
vps[vps<0.2] <- 0.2
vps[vps>0.8] <- 0.8

# sub-units
oms.su <- expand.grid(stk1=names(oms), stk2=names(oms))
nodup <- oms.su
nodup$id <- as.numeric(ac(nodup$stk1))*as.numeric(ac(nodup$stk2))
nodup <- nodup[order(nodup$id),]
nodup <- nodup[nodup$id[-nrow(nodup)]-nodup$id[-1]!=0,]
nodup$id <- paste(nodup$stk1, nodup$stk2, sep=".")

oms.su <- mclapply(split(oms.su, oms.su), function(x){
	stka <- oms[[ac(x$stk1)]]$stk
	stkb <- oms[[ac(x$stk2)]]$stk
	set.seed(seed)
	cres <- rlnorm(prod(dim(catch.n(stka))), 0, log(cvc^2 + 1))
	catch.n(stka) <- catch.n(stka)*cres 
	catch.n(stkb) <- catch.n(stkb)*cres 
	stk <- merge(stka, stkb)
	popspl <- catch.n(stk)
	popspl[] <- vps
	stkc1 <- stk
	stock.n(stkc1) <- popspl*(stock.n(stk))
	catch.n(stkc1) <- catch.n(stk)/stock.n(stk)*stock.n(stkc1)
	catch(stkc1) <- computeCatch(stkc1)
	stkc2 <- stk
	stock.n(stkc2) <- (1-popspl)*(stock.n(stk))
	catch.n(stkc2) <- catch.n(stk)/stock.n(stk)*stock.n(stkc2)
	catch(stkc2) <- computeCatch(stkc2)
	set.seed(seed)
	idxc1 <- getIdx(stkc1, cv=cvq)
	idxc1 <- setPG(idxc1, getPG(stkc1))
	set.seed(seed)
	idxc2 <- getIdx(stkc2, cv=cvq)
	idxc2 <- setPG(idxc2, max(getPG(stkc2), 5))
	list(stkc1=stkc1, idxc1=idxc1, stkc2=stkc2, idxc2=idxc2)
}, mc.cores=mc)

oms.su <- oms.su[nodup$id]

#--------------------------------------------------------------------
# fit
#--------------------------------------------------------------------

fits <- mclapply(oms.su, function(x){
cat("o ") 
	fmod <- ~te(age, year, k = c(10, 10), bs = "tp") + s(year, by=as.numeric(age==1), k=5)
	qmod <- list(~s(age, k=11))
	x$fitc1 <- a4aSCA(x[["stkc1"]], FLIndices(a=x[["idxc1"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x$fitc2 <- a4aSCA(x[["stkc2"]], FLIndices(a=x[["idxc2"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x
}, mc.cores=mc)

# check failures
lst <- lapply(fits, function(x) c(anyNA(harvest(x$fitc1)), anyNA(harvest(x$fitc2))))
fails <- do.call("rbind", lst)

#--------------------------------------------------------------------
# Post-process
#--------------------------------------------------------------------

spatialHeterogeneityRes <- lapply(fits, function(x){
	stks  <- with(x, FLStocks(stkc1=stkc1+fitc1, stkc2=stkc2+fitc2))
	stks$stksmx <- merge(stks$stkc1, stks$stkc2)
	gc()
	list(stks=stks)
})

# recover overall stocks
for(i in 1:length(spatialHeterogeneityRes)) spatialHeterogeneityRes[[i]]$stks$stk <- spatialHeterogeneityRes0[[i]]$stks$stk

diffRes <- mclapply(spatialHeterogeneityRes, function(x){
	ssbmsemx <- with(x, median((ssb(stks$stk)-ssb(stks$stksmx))^2))
	recmsemx <- with(x, median((rec(stks$stk)-rec(stks$stksmx))^2))
	data.frame(ssbmsemx=ssbmsemx, recmsemx=recmsemx)
}, mc.cores=mc)

diffRes <- cbind(do.call("rbind", diffRes), stk1=as.numeric(ac(nodup$stk1)), stk2=as.numeric(ac(nodup$stk2)))

save(spatialHeterogeneityRes, diffRes, vps, oms.su, fits, fails, file="RData.sim2")
rm(spatialHeterogeneityRes, diffRes, oms.su, fits, fails)

#--------------------------------------------------------------------
# diffusion 0.3-0.7
#--------------------------------------------------------------------

# pop split
age <- as.numeric(dimnames(oms[[1]]$stk@catch.n)$age)
vps <- gle(age, b=1, q=6, m=0, nu=0.1)
vps[vps<0.3] <- 0.3
vps[vps>0.7] <- 0.7

# sub-units
oms.su <- expand.grid(stk1=names(oms), stk2=names(oms))
nodup <- oms.su
nodup$id <- as.numeric(ac(nodup$stk1))*as.numeric(ac(nodup$stk2))
nodup <- nodup[order(nodup$id),]
nodup <- nodup[nodup$id[-nrow(nodup)]-nodup$id[-1]!=0,]
nodup$id <- paste(nodup$stk1, nodup$stk2, sep=".")

oms.su <- mclapply(split(oms.su, oms.su), function(x){
	stka <- oms[[ac(x$stk1)]]$stk
	stkb <- oms[[ac(x$stk2)]]$stk
	set.seed(seed)
	cres <- rlnorm(prod(dim(catch.n(stka))), 0, log(cvc^2 + 1))
	catch.n(stka) <- catch.n(stka)*cres 
	catch.n(stkb) <- catch.n(stkb)*cres 
	stk <- merge(stka, stkb)
	popspl <- catch.n(stk)
	popspl[] <- vps
	stkc1 <- stk
	stock.n(stkc1) <- popspl*(stock.n(stk))
	catch.n(stkc1) <- catch.n(stk)/stock.n(stk)*stock.n(stkc1)
	catch(stkc1) <- computeCatch(stkc1)
	stkc2 <- stk
	stock.n(stkc2) <- (1-popspl)*(stock.n(stk))
	catch.n(stkc2) <- catch.n(stk)/stock.n(stk)*stock.n(stkc2)
	catch(stkc2) <- computeCatch(stkc2)
	set.seed(seed)
	idxc1 <- getIdx(stkc1, cv=cvq)
	idxc1 <- setPG(idxc1, getPG(stkc1))
	set.seed(seed)
	idxc2 <- getIdx(stkc2, cv=cvq)
	idxc2 <- setPG(idxc2, max(getPG(stkc2), 5))
	list(stkc1=stkc1, idxc1=idxc1, stkc2=stkc2, idxc2=idxc2)
}, mc.cores=mc)

oms.su <- oms.su[nodup$id]

#--------------------------------------------------------------------
# fit
#--------------------------------------------------------------------

fits <- mclapply(oms.su, function(x){
cat("o ") 
	fmod <- ~te(age, year, k = c(10, 10), bs = "tp") + s(year, by=as.numeric(age==1), k=5)
	qmod <- list(~s(age, k=11))
	x$fitc1 <- a4aSCA(x[["stkc1"]], FLIndices(a=x[["idxc1"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x$fitc2 <- a4aSCA(x[["stkc2"]], FLIndices(a=x[["idxc2"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x
}, mc.cores=mc)

# check failures
lst <- lapply(fits, function(x) c(anyNA(harvest(x$fitc1)), anyNA(harvest(x$fitc2))))
fails <- do.call("rbind", lst)

#--------------------------------------------------------------------
# Post-process
#--------------------------------------------------------------------

spatialHeterogeneityRes <- lapply(fits, function(x){
	stks  <- with(x, FLStocks(stkc1=stkc1+fitc1, stkc2=stkc2+fitc2))
	stks$stksmx <- merge(stks$stkc1, stks$stkc2)
	gc()
	list(stks=stks)
})

# recover overall stocks
for(i in 1:length(spatialHeterogeneityRes)) spatialHeterogeneityRes[[i]]$stks$stk <- spatialHeterogeneityRes0[[i]]$stks$stk

diffRes <- mclapply(spatialHeterogeneityRes, function(x){
	ssbmsemx <- with(x, median((ssb(stks$stk)-ssb(stks$stksmx))^2))
	recmsemx <- with(x, median((rec(stks$stk)-rec(stks$stksmx))^2))
	data.frame(ssbmsemx=ssbmsemx, recmsemx=recmsemx)
}, mc.cores=mc)

diffRes <- cbind(do.call("rbind", diffRes), stk1=as.numeric(ac(nodup$stk1)), stk2=as.numeric(ac(nodup$stk2)))

save(spatialHeterogeneityRes, diffRes, vps, oms.su, fits, fails, file="RData.sim3")
rm(spatialHeterogeneityRes, diffRes, oms.su, fits, fails)

#--------------------------------------------------------------------
# diffusion 0.4-0.6
#--------------------------------------------------------------------

# pop split
age <- as.numeric(dimnames(oms[[1]]$stk@catch.n)$age)
vps <- gle(age, b=1, q=6, m=0, nu=0.1)
vps[vps<0.4] <- 0.4
vps[vps>0.6] <- 0.6

# sub-units
oms.su <- expand.grid(stk1=names(oms), stk2=names(oms))
nodup <- oms.su
nodup$id <- as.numeric(ac(nodup$stk1))*as.numeric(ac(nodup$stk2))
nodup <- nodup[order(nodup$id),]
nodup <- nodup[nodup$id[-nrow(nodup)]-nodup$id[-1]!=0,]
nodup$id <- paste(nodup$stk1, nodup$stk2, sep=".")

oms.su <- mclapply(split(oms.su, oms.su), function(x){
	stka <- oms[[ac(x$stk1)]]$stk
	stkb <- oms[[ac(x$stk2)]]$stk
	set.seed(seed)
	cres <- rlnorm(prod(dim(catch.n(stka))), 0, log(cvc^2 + 1))
	catch.n(stka) <- catch.n(stka)*cres 
	catch.n(stkb) <- catch.n(stkb)*cres 
	stk <- merge(stka, stkb)
	popspl <- catch.n(stk)
	popspl[] <- vps
	stkc1 <- stk
	stock.n(stkc1) <- popspl*(stock.n(stk))
	catch.n(stkc1) <- catch.n(stk)/stock.n(stk)*stock.n(stkc1)
	catch(stkc1) <- computeCatch(stkc1)
	stkc2 <- stk
	stock.n(stkc2) <- (1-popspl)*(stock.n(stk))
	catch.n(stkc2) <- catch.n(stk)/stock.n(stk)*stock.n(stkc2)
	catch(stkc2) <- computeCatch(stkc2)
	set.seed(seed)
	idxc1 <- getIdx(stkc1, cv=cvq)
	idxc1 <- setPG(idxc1, getPG(stkc1))
	set.seed(seed)
	idxc2 <- getIdx(stkc2, cv=cvq)
	idxc2 <- setPG(idxc2, max(getPG(stkc2), 5))
	list(stkc1=stkc1, idxc1=idxc1, stkc2=stkc2, idxc2=idxc2)
}, mc.cores=mc)

oms.su <- oms.su[nodup$id]

#--------------------------------------------------------------------
# fit
#--------------------------------------------------------------------

fits <- mclapply(oms.su, function(x){
cat("o ") 
	fmod <- ~te(age, year, k = c(10, 10), bs = "tp") + s(year, by=as.numeric(age==1), k=5)
	qmod <- list(~s(age, k=11))
	x$fitc1 <- a4aSCA(x[["stkc1"]], FLIndices(a=x[["idxc1"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x$fitc2 <- a4aSCA(x[["stkc2"]], FLIndices(a=x[["idxc2"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x
}, mc.cores=mc)

# check failures
lst <- lapply(fits, function(x) c(anyNA(harvest(x$fitc1)), anyNA(harvest(x$fitc2))))
fails <- do.call("rbind", lst)

#--------------------------------------------------------------------
# Post-process
#--------------------------------------------------------------------

spatialHeterogeneityRes <- lapply(fits, function(x){
	stks  <- with(x, FLStocks(stkc1=stkc1+fitc1, stkc2=stkc2+fitc2))
	stks$stksmx <- merge(stks$stkc1, stks$stkc2)
	gc()
	list(stks=stks)
})

# recover overall stocks
for(i in 1:length(spatialHeterogeneityRes)) spatialHeterogeneityRes[[i]]$stks$stk <- spatialHeterogeneityRes0[[i]]$stks$stk

diffRes <- mclapply(spatialHeterogeneityRes, function(x){
	ssbmsemx <- with(x, median((ssb(stks$stk)-ssb(stks$stksmx))^2))
	recmsemx <- with(x, median((rec(stks$stk)-rec(stks$stksmx))^2))
	data.frame(ssbmsemx=ssbmsemx, recmsemx=recmsemx)
}, mc.cores=mc)

diffRes <- cbind(do.call("rbind", diffRes), stk1=as.numeric(ac(nodup$stk1)), stk2=as.numeric(ac(nodup$stk2)))

save(spatialHeterogeneityRes, diffRes, vps, oms.su, fits, fails, file="RData.sim4")
rm(spatialHeterogeneityRes, diffRes, oms.su, fits, fails)

#--------------------------------------------------------------------
# diffusion 1-0
#--------------------------------------------------------------------

# pop split
age <- as.numeric(dimnames(oms[[1]]$stk@catch.n)$age)
vps <- gle(age, b=1, q=6, m=0, nu=0.1)

# sub-units
oms.su <- expand.grid(stk1=names(oms), stk2=names(oms))
nodup <- oms.su
nodup$id <- as.numeric(ac(nodup$stk1))*as.numeric(ac(nodup$stk2))
nodup <- nodup[order(nodup$id),]
nodup <- nodup[nodup$id[-nrow(nodup)]-nodup$id[-1]!=0,]
nodup$id <- paste(nodup$stk1, nodup$stk2, sep=".")

oms.su <- mclapply(split(oms.su, oms.su), function(x){
	stka <- oms[[ac(x$stk1)]]$stk
	stkb <- oms[[ac(x$stk2)]]$stk
	set.seed(seed)
	cres <- rlnorm(prod(dim(catch.n(stka))), 0, log(cvc^2 + 1))
	catch.n(stka) <- catch.n(stka)*cres 
	catch.n(stkb) <- catch.n(stkb)*cres 
	stk <- merge(stka, stkb)
	popspl <- catch.n(stk)
	popspl[] <- vps
	stkc1 <- stk
	stock.n(stkc1) <- popspl*(stock.n(stk))
	catch.n(stkc1) <- catch.n(stk)/stock.n(stk)*stock.n(stkc1)
	catch(stkc1) <- computeCatch(stkc1)
	stkc2 <- stk
	stock.n(stkc2) <- (1-popspl)*(stock.n(stk))
	catch.n(stkc2) <- catch.n(stk)/stock.n(stk)*stock.n(stkc2)
	catch(stkc2) <- computeCatch(stkc2)
	set.seed(seed)
	idxc1 <- getIdx(stkc1, cv=cvq)
	idxc1 <- setPG(idxc1, getPG(stkc1))
	set.seed(seed)
	idxc2 <- getIdx(stkc2, cv=cvq)
	idxc2 <- setPG(idxc2, max(getPG(stkc2), 5))
	list(stkc1=stkc1, idxc1=idxc1, stkc2=stkc2, idxc2=idxc2)
}, mc.cores=mc)

oms.su <- oms.su[nodup$id]

#--------------------------------------------------------------------
# fit
#--------------------------------------------------------------------

fits <- mclapply(oms.su, function(x){
cat("o ") 
	fmod <- ~te(age, year, k = c(10, 10), bs = "tp") + s(year, by=as.numeric(age==1), k=5)
	qmod <- list(~s(age, k=11))
	x$fitc1 <- a4aSCA(x[["stkc1"]], FLIndices(a=x[["idxc1"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x$fitc2 <- a4aSCA(x[["stkc2"]], FLIndices(a=x[["idxc2"]]), fit="MP", fmodel=fmod , qmodel=qmod)#, center=1)
	x
}, mc.cores=mc)

# check failures
lst <- lapply(fits, function(x) c(anyNA(harvest(x$fitc1)), anyNA(harvest(x$fitc2))))
fails <- do.call("rbind", lst)

#--------------------------------------------------------------------
# Post-process
#--------------------------------------------------------------------

spatialHeterogeneityRes <- lapply(fits, function(x){
	stks  <- with(x, FLStocks(stkc1=stkc1+fitc1, stkc2=stkc2+fitc2))
	stks$stksmx <- merge(stks$stkc1, stks$stkc2)
	gc()
	list(stks=stks)
})

# recover overall stocks
for(i in 1:length(spatialHeterogeneityRes)) spatialHeterogeneityRes[[i]]$stks$stk <- spatialHeterogeneityRes0[[i]]$stks$stk

diffRes <- mclapply(spatialHeterogeneityRes, function(x){
	ssbmsemx <- with(x, median((ssb(stks$stk)-ssb(stks$stksmx))^2))
	recmsemx <- with(x, median((rec(stks$stk)-rec(stks$stksmx))^2))
	data.frame(ssbmsemx=ssbmsemx, recmsemx=recmsemx)
}, mc.cores=mc)

diffRes <- cbind(do.call("rbind", diffRes), stk1=as.numeric(ac(nodup$stk1)), stk2=as.numeric(ac(nodup$stk2)))

save(spatialHeterogeneityRes, diffRes, vps, oms.su, fits, fails, file="RData.sim5")
rm(spatialHeterogeneityRes, diffRes, oms.su, fits, fails)

#====================================================================
# putting all together
#====================================================================
attach("RData.sim0")
diffRes$diffusion <- 0
diffResAll <- diffRes
rm(diffRes)
detach()
attach("RData.sim1")
diffRes$diffusion <- 0.8
diffRes$ssbmseu <- diffRes$recmseu <- NA
diffResAll <- rbind(diffResAll, diffRes[,names(diffResAll)])
rm(diffRes)
detach()
attach("RData.sim2")
diffRes$diffusion <- 0.6
diffRes$ssbmseu <- diffRes$recmseu <- NA
diffResAll <- rbind(diffResAll, diffRes[,names(diffResAll)])
detach()
rm(diffRes)
attach("RData.sim3")
diffRes$diffusion <- 0.4
diffRes$ssbmseu <- diffRes$recmseu <- NA
diffResAll <- rbind(diffResAll, diffRes[,names(diffResAll)])
detach()
rm(diffRes)
attach("RData.sim4")
diffRes$diffusion <- 0.2
diffRes$ssbmseu <- diffRes$recmseu <- NA
diffResAll <- rbind(diffResAll, diffRes[,names(diffResAll)])
rm(diffRes)
detach()
attach("RData.sim5")
diffRes$diffusion <- 1
diffRes$ssbmseu <- diffRes$recmseu <- NA
diffResAll <- rbind(diffResAll, diffRes[,names(diffResAll)])
rm(diffRes)
detach()

save(diffResAll, file="RData.simsumm")
#--------------------------------------------------------------------
# plotting
#--------------------------------------------------------------------
load("RData.simsumm")
df0 <- melt(diffResAll, id.vars=5:7)
v <- as.character(df0$variable)
df0$sce <- substr(v, nchar(v), nchar(v))
df0$stat <- substr(v, 1, 3)
df0 <- transform(df0, phat= 0.5*(stk1+stk2), rat= stk1/stk2, diff= abs(stk1-stk2), sce=factor(sce, labels=c("I","D")))
df0 <- subset(df0, rat < 1.01)
df0$stat <- factor(df0$stat, labels=c("R","SSB"))
df0$grp <- paste(df0$diffusion, df0$sce, sep=":") 
df0 <- df0[!is.na(df0$value),]

pdf("msessbd.pdf", 10, 10)
xyplot(sqrt(value)~rat, groups=diffusion, type=c("smooth", "p"), auto.key=list(columns=1, space="right"), ylab="mean square error", xlab="steepness ratio", data=subset(df0, stat=="SSB" & sce=="D"), scales=list(y=list(relation="free")), pch="o", main="Comparison of sqrt(MSE) of SSB between disagregated and aggregated assessments for distinct diffusion processes. \n A value of 0.5 (in green) refers to a weak diffusion where the population was split in half. \n While a value of 0 (in blue) refers to the strongest \n diffusion process, where all recruits are in one area and all adults in another area.")
dev.off()

pdf("msessbdVSi.pdf", 10, 10)
xyplot(sqrt(value)~rat, groups=sce, type=c("smooth", "p"), auto.key=list(columns=1, space="right"), ylab="mean square error", xlab="steepness ratio", data=subset(df0, stat=="SSB" & grp %in% c("0:I", "0:D") ), scales=list(y=list(relation="free")), pch="o", main="Comparison of sqrt(MSE) of SSB between disagregated and aggregated assessments for distinct scenarios of non-diffusion. \n Scenario I is build of merging two independent stocks, although with different steepness (in the x axis). \n Scenario D is build by using a flat diffusion (green line in logistic plot), which generates sub-units with similar dinamics.")
dev.off()

# pop split
age <- vps <- 0:maxage
vps[] <- 1
vps[1] <- 0
df0 <- data.frame(age=age, vps=vps, diffusion=1)

df1 <- df0 
df1$vps[df1$vps<0.1] <- 0.1
df1$vps[df1$vps>0.9] <- 0.9
df1$diffusion <- 0.9-0.1
df0 <- rbind(df0, df1)

df1$vps[df1$vps<0.2] <- 0.2
df1$vps[df1$vps>0.8] <- 0.8
df1$diffusion <- 0.8-0.2
df0 <- rbind(df0, df1)

df1$vps[df1$vps<0.3] <- 0.3
df1$vps[df1$vps>0.7] <- 0.7
df1$diffusion <- 0.7-0.3
df0 <- rbind(df0, df1)

df1$vps[df1$vps<0.4] <- 0.4
df1$vps[df1$vps>0.6] <- 0.6
df1$diffusion <- 0.6-0.4
df0 <- rbind(df0, df1)

df1$vps <- 0.5
df1$diffusion <- 0
diffusionModel <- rbind(df0, df1)

pdf("diffproc.pdf", 10, 10)
xyplot(vps~age, groups=diffusion, data=diffusionModel, type="l", auto.key=list(points=FALSE, lines=TRUE, space="right"), main="Logistic models used for simulating the diffusion process. \n A value of 0.5 generates a flat line (in green) which means the population was split in half. \n While a value of 0 generates a logistic between 0 and 1 (in blue) which imposes the strongest \n diffusion process, where all recruits are in one area and all adults in another area.")
dev.off()





