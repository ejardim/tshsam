###############################################################################
# EJ(20120413)
# Auxiliary functions for LH generation of datasets
# Parameters: 
#	a1, sr, sl = 50% selectivity age, right variance, left variance
#	s, v = steepness and virgin biomass for S/R
#	M1, M2 = two components of M
#	a50, asym = age of 50% maturity, ???
#	linf, k, t0 = vonBertalanffy growth pars
#	a, b = length~weight pars
#	bg ?? ato95 ??
###############################################################################

iterMedians <- function(x, ...){
	return(apply(x, c(1:5), median, na.rm = FALSE))
}

iterSums <- function(x, ...){
	return(apply(x, c(1:5), sum, na.rm = FALSE))
}

iterCv <- function(object, ...){
	sqrt(iterVars(object))/iterMeans(object)
}

#==============================================================================
# gislasim - cleaned version
#==============================================================================

setGeneric("gislasim", function(linf, ...) standardGeneric("gislasim"))

setMethod("gislasim", signature(linf="numeric"), function (linf, t0 = -0.1, a = 1e-05, b = 3, ato95 = 1, sl = 2, sr = 5000, s = 0.9, v = 1000, asym=1, bg=b, iter=1, k="missing", M1="missing", M2="missing", a50="missing", a1="missing"){
    if(missing(k))  k <- 3.15 * linf^(-0.64)
    if(missing(M1)) M1 <- 0.55 + 1.44 * log(linf) + log(k) 
    if(missing(M2)) M2 <- -1.61
    if(missing(a50)) a50 <- FLBRP:::invVonB(FLPar(linf=linf, t0=t0, k=k), 0.72 * linf^0.93)
    if(missing(a1)) a1 <- a50
    par <- FLPar(linf=linf, k=k, t0 = t0, a = a, b = b, asym=asym, bg=bg, sl=sl, sr=sr, s=s, v=v, M1=M1, M2=M2, ato95 = ato95, a50=a50, a1=a1, iter=iter)
    attributes(par)$units = c("cm", "kg", "1000s")
    return(par)
})

setMethod("gislasim", signature(linf="FLPar"), function (linf){
    # Renaming to avoid confusing the argument with the object.
    # linf here is an FLPar object that can contain several parameters 
    object <- linf
    rm(linf)
    # now the real thing
    v0 <- dimnames(object)$params	    
    if(!("linf" %in% v0)) stop("The function requires linf.")
    par <- FLPar(c(linf=NA, t0 = -0.1, a = 1e-05, b = 3, ato95 = 1, sl = 2, sr = 5000, s = 0.9, v = 1000, asym=1, bg=3, k=NA, M1=NA, M2=NA, a50=NA, a1=NA), iter=ncol(object))
    dimnames(par)$iter <- dimnames(object)$iter 
    par[dimnames(object)$params] <- object
    if(!("bg" %in% v0)) par["bg"] = par["b"]
    if(!("k" %in% v0)) par["k"] = 3.15 * par["linf"]^(-0.64)
    if(!("M1" %in% v0)) par["M1"] = 0.55 + 1.44 * log(par["linf"]) + log(par["k"])
    if(!("M2" %in% v0)) par["M2"] = -1.61
    if(!("a50" %in% v0)) par["a50"] = FLBRP:::invVonB(FLPar(linf=par["linf"], t0=par["t0"], k=par["k"]), c(0.72 * par["linf"]^0.93))
    if(!("a1" %in% v0)) par["a1"] = par["a50"]
    attributes(par)$units = c("cm", "kg", "1000s")
    return(par)
})

#==============================================================================
# genFunctions
#==============================================================================

genObs <- function(stock, qcv=0.1, ccv=0.1) {
	ages <- 1:range(stock)["maxfbar"]
	n <- stock.n(stock)[ac(ages)]
	z <- harvest(stock)[ac(ages)] + m(stock)[ac(ages)]
	logq <- -exp(-exp(0.2 * ages)) - 3 # trawl like catchability
	# observe index in 1st quarter with qcv
	index <- FLIndex(index = n[ac(ages)] * exp(-0.25 * z[ac(ages)]) * exp(logq + rnorm(prod(dim(n[ac(ages)])), 0, qcv)))  # 10% cv
	# OR observe index in 1st quarter with qcv and cov N
	#Sig2 <- cor(t(n[drop = T]))*qcv^2
	#index <- FLIndex(index = n * exp(-0.25 * z) * c(exp(logq + t(mvrnorm(dim(n)[2], rep(0,dim(n)[1]), Sig2)))))
	range(index)[c("startf","endf")] <- 0.25
	# observe catch with ccv
	catch.n(stock) <- catch.n(stock) * exp(rnorm(prod(dim(catch.n(stock))), 0, ccv))
	catch(stock) <- computeCatch(stock)
	list(stock = stock, index = list(index))
}

genIdx <- function(stock, qcv=0.1, ccv=0.1) {
	ages <- 1:range(stock)["maxfbar"]
	n <- stock.n(stock)[ac(ages)]
	z <- harvest(stock)[ac(ages)] + m(stock)[ac(ages)]
	logq <- -exp(-exp(0.2 * ages)) - 3 # trawl like catchability
	# observe index in 1st quarter with qcv
	index <- FLIndex(index = n[ac(ages)] * exp(-0.25 * z[ac(ages)]) * exp(logq + rnorm(prod(dim(n[ac(ages)])), 0, qcv)))  # 10% cv
	# OR observe index in 1st quarter with qcv and cov N
	#Sig2 <- cor(t(n[drop = T]))*qcv^2
	#index <- FLIndex(index = n * exp(-0.25 * z) * c(exp(logq + t(mvrnorm(dim(n)[2], rep(0,dim(n)[1]), Sig2)))))
	range(index)[c("startf","endf")] <- 0.25
	list(stock = stock, index = list(index))
}

#==============================================================================
# create index
#==============================================================================

setGeneric("getIdx", function(object, ...) standardGeneric("getIdx"))

setMethod("getIdx", signature(object="FLStock"), function (object, q=0.01, cv=0.1, ...){
	flq <- object@stock.n	
	age <- as.numeric(dimnames(flq)$age)
	sel <- exp(-age/2)
	sel <- sel/max(sel)
	flq[] <- sel%o%(q*rlnorm(length(dimnames(flq)$year), 0, log(cv^2+1)))
	idx <- FLIndex(index=stock.n(object)*flq)	
	range(idx)[c("startf", "endf")] <- c(0, 0)
	idx
	
})

#==============================================================================
# merge stocks
#==============================================================================

setGeneric("merge", function(x, y, ...) standardGeneric("merge"))

setMethod("merge", signature(x="FLStock", y="FLStock"), function (x, y, ...){

	# catch.n
	lst <- mcf(list(x@catch.n, y@catch.n))
	lst <- lapply(lst, function(x){ x[is.na(x)] <- 0; x})
	stk <- FLStock(catch.n=lst[[1]])
	stk@catch.n <- lst[[1]] + lst[[2]]
	lstwt <- mcf(list(x@catch.wt, y@catch.wt))
	lstwt <- lapply(lstwt, function(x){ x[is.na(x)] <- max(x, na.rm=TRUE); x})
	stk@catch.wt <- 1/stk@catch.n * (lstwt[[1]]*lst[[1]] + lstwt[[2]]*lst[[2]])
	stk@catch <- computeCatch(stk)
	
	# lstwt <- mcf(list(x@harvest.spwn, y@harvest.spwn))
	# lstwt <- lapply(lstwt, function(x){ x[is.na(x)] <- 0; x})
	# stk@harvest.spwn <- 1/stk@catch.n * (lstwt[[1]]*lst[[1]] + lstwt[[2]]*lst[[2]])
	stk@harvest.spwn[] <- 0
	
	# stock.n
	lst <- mcf(list(x@stock.n, y@stock.n))
	lst <- lapply(lst, function(x){ x[is.na(x)] <- 0; x})
	stk@stock.n <- lst[[1]] + lst[[2]]
	lstwt <- mcf(list(x@stock.wt, y@stock.wt))
	lstwt <- lapply(lstwt, function(x){ x[is.na(x)] <- max(x, na.rm=TRUE); x})
	stk@stock.wt <- 1/stk@stock.n * (lstwt[[1]]*lst[[1]] + lstwt[[2]]*lst[[2]])
	stk@stock <- computeStock(stk)

	# mat
	lstwt <- mcf(list(x@mat, y@mat))
	lstwt <- lapply(lstwt, function(x){ x[is.na(x)] <- max(x, na.rm=TRUE); x})
	stk@mat <- 1/stk@stock.n * (lstwt[[1]]*lst[[1]] + lstwt[[2]]*lst[[2]])

	# m
	lstwt <- mcf(list(x@m, y@m))
	lstwt <- lapply(lstwt, function(x){ x[is.na(x)] <- max(x, na.rm=TRUE); x})
	stk@m <- 1/stk@stock.n * (lstwt[[1]]*lst[[1]] + lstwt[[2]]*lst[[2]])

	# f
	harvest(stk)[-dims(stk)$age, -dims(stk)$year] <- -log(stock.n(stk)[-1,-1]/stock.n(stk)[-dims(stk)$age, -dims(stk)$year]) - m(stk)[-dims(stk)$age, -dims(stk)$year]
	harvest(stk)[dims(stk)$age] <- harvest(stk)[dims(stk)$age-1] <- harvest(stk)[dims(stk)$age-2]
	harvest(stk)[, dims(stk)$year] <- harvest(stk)[, dims(stk)$year-1]
	units(stk@harvest) <- units(x@harvest)

#	lstwt <- mcf(list(x@m.spwn, y@m.spwn))
#	lstwt <- lapply(lstwt, function(x){ x[is.na(x)] <- 0; x})
#	stk@m.spwn <- 1/stk@stock.n * (lstwt[[1]]*lst[[1]] + lstwt[[2]]*lst[[2]])
	stk@m.spwn[] <- 0
	
	# landings
	lst <- mcf(list(x@landings.n, y@landings.n))
	lst <- lapply(lst, function(x){ x[is.na(x)] <- 0; x})
	stk@landings.n <- lst[[1]] + lst[[2]]
	lstwt <- mcf(list(x@landings.wt, y@landings.wt))
	lstwt <- lapply(lstwt, function(x){ x[is.na(x)] <- max(x, na.rm=TRUE); x})
	stk@landings.wt <- 1/stk@landings.n * (lstwt[[1]]*lst[[1]] + lstwt[[2]]*lst[[2]])
	stk@landings <- computeLandings(stk)

	# discards
	lst <- mcf(list(x@discards.n, y@discards.n))
	lst <- lapply(lst, function(x){ x[is.na(x)] <- 0; x})
	stk@discards.n <- lst[[1]] + lst[[2]]
	lstwt <- mcf(list(x@discards.wt, y@discards.wt))
	lstwt <- lapply(lstwt, function(x){ x[is.na(x)] <- max(x, na.rm=TRUE); x})
	stk@discards.wt <- 1/stk@discards.n * (lstwt[[1]]*lst[[1]] + lstwt[[2]]*lst[[2]])
	stk@discards <- computeDiscards(stk)

	stk	

})

#==============================================================================
# plus group	
#==============================================================================

setGeneric("setPG", function(object, ...) standardGeneric("setPG"))

setMethod("setPG", signature(object="FLStock"), function (object, threshold=0.99, minPG=6, ...){
	flq <- catch.n(object)
	av <- as.numeric(dimnames(flq)$age)
	flq <- apply(flq, 2, function(x) cumsum(x)/sum(x, na.rm=T))
	setPlusGroup(object, plusgroup=max(min(av[apply(flq, 1, max, na.rm=T)>threshold]), minPG))
})

setMethod("setPG", signature(object="FLIndex"), function (object, plusgroup, ...){
	object <- object[ac(range(object)["min"]:min(range(object)["max"], plusgroup-1))]
	range(object)["plusgroup"] <- NA
	object
})

setGeneric("getPG", function(object, ...) standardGeneric("getPG"))

setMethod("getPG", signature(object="FLStock"), function (object){
	range(object)["plusgroup"]
})

setMethod("getPG", signature(object="FLIndex"), function (object){
	range(object)["plusgroup"]
})

#==============================================================================
# baranov	
#==============================================================================

setGeneric("baranov", function(object, ...) standardGeneric("baranov"))

setMethod("baranov", signature(object="FLStock"), function (object, ...){
	N <- stock.n(object)
	F <- harvest(object)
	M <- m(object)
	F/(F+M)*(1-exp(-(F+M)))*N
})

#==============================================================================
# baranov	
#==============================================================================

setGeneric("getFitStats", function(object, ...) standardGeneric("getFitStats"))

setMethod("getFitStats", signature(object="a4aFit"), function (object, ...){
	data.frame(t(fitSumm(object)), AIC=AIC(object), BIC=BIC(object))
})


#==============================================================================
# ypr	
#==============================================================================

setMethod("ypr", signature(object="FLStock"), function (object, ...){
	catch(object)/rec(object)
})

#==============================================================================
# run assessment and update	
#==============================================================================

runsa <- function(x){
	stk <- x$stk
	dms <- lapply(dimnames(catch.n(stk)), length)
	ka <- dms$age
	ka <- if (ka < 3) ka else min(max(4, ceiling(.5 * ka)), 10)
	ky <- floor(.75 * dms$year)
	if (ka >= 4) {
		fmodel <- formula(paste("~ te(age, year, k = c(", ka,",", ky,"), bs = 'tp')"))
	} else {
		fmodel <- formula(paste("~ factor(age) + s(year, k = ", ky,")"))
	}

	idx <- FLIndices(idx=x$idx)
	fit <- sca(stk, idx, fit="assessment", fmodel=fmodel)
	return(fit)
}

#==============================================================================
# Generalized Logistic	
#==============================================================================

# The generalised logistic function or curve, also known as Richards' curve, 
# originally developed for growth modelling, is an extension of the logistic or 
# sigmoid functions, allowing for more flexible S-shaped curves:
# It has seven parameters: 
#	A: the lower asymptote;
#	K: the upper asymptote. If A=0 then K is called the carrying capacity;
#	B: the growth rate;
#	nu > 0 : affects near which asymptote maximum growth occurs.
#	Q: is related to the value Y(0)
#	M: starting time, t_0 
#	C: typically takes a value of 1.

gle <- function(x, a=0, k=1, b, nu=0.5, q=0, m=0, c=1) a + ((k-a)/(c + q * exp(-b *(x - m)))^(1/nu))


#==============================================================================
# extract	
#==============================================================================

getStkInfo <- function(object, scn, hcr){
	res <- lapply(object, function(x){
		lh <- as.data.frame(attr(x, "lhPars"))
		rp <- as.data.frame(attr(x, "refpts")[c("msy", "f0.1", "fmax"), "harvest"])
		data.frame(params=c(ac(lh[,1]), ac(rp[,1])), value=c(lh[,3], rp[,4]), scn=scn, hcr=hcr)	
	})
	
	res <- do.call("rbind", res)
	res$stk <- unlist(lapply(strsplit(rownames(res), "[.]"), "[[", 1))
	rownames(res) <- NULL
	res
	
}
	
getRes <- function(object, scn, hcr){
	res <- lapply(object, as.data.frame)
	res <- do.call("rbind", res)
	res$stk <- unlist(lapply(strsplit(rownames(res), "[.]"), "[[", 1))
	rownames(res) <- NULL
	res$scn <- scn
	res$hcr <- hcr
	res
}

#==============================================================================
# plot for presentations	
#==============================================================================

plot4presentations <- function (x, main = "", xlab = "", ylab = "", na.rm = TRUE, ...) {
        dup <- duplicated(names(x))
        if (any(dup)) {
            names(x)[dup] <- paste(names(x)[dup], LETTERS[seq(sum(dup))], 
                sep = "_")
            warning("Duplicated names in object, changed to differentiate")
        }
        fqs <- lapply(x, function(y) FLQuants(Rec = rec(y), SSB = ssb(y)))
#        fqs <- lapply(x, function(y) FLQuants(Rec = rec(y), SSB = ssb(y), 
#            Catch = catch(y), Harvest = fbar(y)))
        its <- unlist(lapply(x, function(x) dims(x)$iter))
        if (any(its > 1)) {
            fqs <- lapply(fqs, function(y) as.data.frame(lapply(y, 
                quantile, c(0.33, 0.5, 0.66), na.rm = TRUE)))
        }
        else {
            fqs <- lapply(fqs, as.data.frame)
            fqs <- lapply(fqs, function(x) {
                x$iter <- "50%"
                return(x)
            })
        }
        stk <- rep.int(names(fqs), unlist(lapply(fqs, nrow)))
        fqs <- do.call(rbind, fqs)
        rownames(fqs) <- NULL
        fqs <- transform(fqs, stock = stk)
        df <- dcast(fqs, age + year + unit + season + area + 
            qname + stock ~ iter, value.var = "data")
        p <- ggplot(data = df, aes(x = year, y = `50%`, group = stock)) + 
            facet_grid(. ~ qname, scales = "free") + geom_line(aes(colour = stock), na.rm = na.rm) + xlab(xlab) + ylab(ylab) + expand_limits(y = 0) + theme(legend.title = element_blank(), legend.text=element_text(size=22), panel.background = element_blank(), panel.grid.major=element_blank(), strip.text.x = element_text(size = 18), strip.background=element_rect(fill = 'white'))
        if (any(unlist(lapply(x, function(y) dims(y)$iter)) > 
            1)) {
            p <- p + geom_ribbon(aes(x = year, ymin = `33%`, 
                ymax = `66%`, group = stock, colour = stock, 
                fill = stock), alpha = 0.5, linetype = 0, na.rm = na.rm)
        }
        return(p)
    }

#==============================================================================
# cohort cross validation
#==============================================================================
gcv <- function(stk, ids, cohorts=range(stk)[c("minyear")]:range(stk)[c("maxyear")]){
	coh <- as.character(cohorts)
	cv <- vector(length=length(coh))
	names(cv) <- coh
	cth <- catch.n(stk)
	for(i in coh){
		cat(i, ",")
		flc <- FLCohort(catch.n(stk))
		flc[,i] <- NA
		catch.n(stk) <- as(flc, "FLQuant")
		fit <- sca(stk, ids)
		cv[i] <- mean(log(FLCohort(cth)[,i]/FLCohort(catch.n(fit))[,i])^2, na.rm=T)
		catch.n(stk) <- cth
	}
	cv
}

#==============================================================================
# check F@age shape
#==============================================================================

checkFshape <- function(stock, indices, file="checkf.pdf", ...){
	qmod <- list(~s(age, k = 5), ~s(age, k = 4))
	pdf(file)
	na <- 3:dims(stock)$age
	for(i in na){
		fmod <- as.formula(paste("~ te(age, year, k = c(", i, ", 5))"))
		fit <- sca(stk, ids, fit="assessment", fmodel=fmod, qmodel=qmod)
		res <- residuals(fit, stock, indices)
		print(plot(res))
		print(wireframe(data~year+age, data=harvest(fit), main=i))
	}
	dev.off()
}

#==============================================================================
# permutation
#==============================================================================

setGeneric("quantPerm", function(object, ...) standardGeneric("quantPerm"))

setMethod("quantPerm", signature(object="FLQuant"), function (object, seed="missing", ...){
	if(missing(seed)){
		object[] <- t(apply(object, 1, sample, size=ncol(object), replace=FALSE))
	} else {
		set.seed(seed)
		object[] <- t(apply(object, 1, sample, size=ncol(object), replace=FALSE))
	}
	object
})

#==============================================================================
# permutation test for sca
#==============================================================================

setGeneric("perm", function(stock, indices, ...) standardGeneric("perm"))

setMethod("perm", signature("FLStock", "FLIndices"), function(stock, indices, fmodel, qmodel, srmodel = ~ factor(year), fit = "assessment", n=10, seed="missing"){
	flq <- catch.n(stock)
	flq <- flq[,,,,,rep(1, n)]
	if(missing(seed)){
		flq[,,1,1,1,] <- aperm(apply(flq, c(1, 6), sample, size=ncol(flq), replace=FALSE), c(2,1,3))
	} else {
		set.seed(seed)
		flq[,,1,1,1,] <- aperm(apply(flq, c(1, 6), sample, size=ncol(flq), replace=FALSE), c(2,1,3))
	}
	flq -> catch.n(stock)
	stock + sca(stock, indices, fmodel, qmodel, fit=fit)
})




