estimateCellCounts <- function (rgSet, compositeCellType = "Blood", 
								processMethod = "auto", probeSelect = "auto",
                                cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"),
                                returnAll = FALSE, meanPlot = FALSE, verbose=TRUE, ...) {
    platform <- annotation(rgSet)[which(names(annotation(rgSet))=="array")]
	platform <- if(platform == "IlluminaHumanMethylation450k") {platform <- "450k"}
	if((compositeCellType == "CordBlood") && (!"nRBC" %in% cellTypes))
		cat("[estimateCellCounts] Consider including 'nRBC' in argument 'cellTypes' for cord blood estimation.\n")   
	referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
    subverbose <- max(as.integer(verbose) - 1L, 0L)
    if(!require(referencePkg, character.only = TRUE))
        stop(sprintf("Could not find reference data package for compositeCellType '%s' and platform '%s' (inferred package name is '%s')",
                     compositeCellType, platform, referencePkg))
    data(list = referencePkg) 
    referenceRGset <- get(referencePkg) 
    if(! "CellType" %in% names(pData(referenceRGset)))
        stop(sprintf("the reference sorted dataset (in this case '%s') needs to have a phenoData column called 'CellType'"),
             names(referencePkg))
    if(sum(colnames(rgSet) %in% colnames(referenceRGset)) > 0)
        stop("the sample/column names in the user set must not be in the reference data ")
    if(!all(cellTypes %in% referenceRGset$CellType))
        stop(sprintf("all elements of argument 'cellTypes' needs to be part of the reference phenoData columns 'CellType' (containg the following elements: '%s')",
                     paste(unique(referenceRGset$cellType), collapse = "', '")))
	if(length(unique(cellTypes)) < 2)
        stop("At least 2 cell types must be provided.")
	if ((processMethod == "auto") && (compositeCellType == "CordBlood")){
		processMethod <- "Noob"} 
	if ((processMethod == "auto") && (compositeCellType != "CordBlood")){
		processMethod <- "Quantile"}
	processMethod <- get(list(sprintf("preprocess%s",processMethod))[[1]])
	if ((probeSelect == "auto") && (compositeCellType == "CordBlood")){
		probeSelect <- "any"} 
	if ((probeSelect == "auto") && (compositeCellType != "CordBlood")){
		probeSelect <- "both"} 	

    if(verbose) cat("[estimateCellCounts] Combining user data with reference (flow sorted) data.\n")
    newpd <- data.frame(sampleNames = c(sampleNames(rgSet), sampleNames(referenceRGset)),
                        studyIndex = rep(c("user", "reference"),
                        times = c(ncol(rgSet), ncol(referenceRGset))),
                        stringsAsFactors=FALSE)
	referencePd <- pData(referenceRGset)
	pData(referenceRGset) <- data.frame() #To avoid errors in the combine call
	combinedRGset <- combine(rgSet, referenceRGset)
    pData(combinedRGset) <- newpd
	sampleNames(combinedRGset) <- newpd$sampleNames #Avoid errors in .isRG call in Funnorm
    rm(referenceRGset)
    
	if(verbose) cat("[estimateCellCounts] Processing user and reference data together.\n")
	if (compositeCellType=="CordBlood"){  
    combinedMset <- processMethod(combinedRGset,verbose=subverbose)
	compTable <- get(paste(referencePkg,".compTable",sep="")) 
	combinedMset <- combinedMset[which(rownames(combinedMset)%in%rownames(compTable)),]
    rm(combinedRGset)} else {
    combinedMset <- processMethod(combinedRGset) 
    rm(combinedRGset)}
 
    ## Extracts normalized reference data 
    referenceMset <- combinedMset[, combinedMset$studyIndex == "reference"]
    pData(referenceMset) <- as(referencePd, "DataFrame")
    mSet <- combinedMset[, combinedMset$studyIndex == "user"]
	pData(mSet) <- as(pData(rgSet), "DataFrame")
    rm(combinedMset)
    
    if(verbose) cat("[estimateCellCounts] Picking probes for composition estimation.\n")
    compData <- pickCompProbes(referenceMset, cellTypes = cellTypes, compositeCellType = compositeCellType, probeSelect = probeSelect)
    coefs <- compData$coefEsts
    rm(referenceMset)
    
    if(verbose) cat("[estimateCellCounts] Estimating composition.\n")
    counts <- projectCellType(getBeta(mSet)[rownames(coefs), ], coefs)
    rownames(counts) <- sampleNames(rgSet)

    if (meanPlot) {
        smeans <- compData$sampleMeans
        smeans <- smeans[order(names(smeans))]
        sampleMeans <- c(colMeans(getBeta(mSet)[rownames(coefs), ]), smeans)

        sampleColors <- c(rep(1, ncol(mSet)), 1 + as.numeric(factor(names(smeans))))
        plot(sampleMeans, pch = 21, bg = sampleColors)
        legend("bottomleft", c("blood", levels(factor(names(smeans)))),
               col = 1:7, pch = 15)
    }
    if(returnAll) {
        list(counts = counts, compTable = compData$compTable,
             normalizedData = mSet)
    } else {
        counts
    }
}

pickCompProbes <- function(mSet, cellTypes = NULL, numProbes = 50, compositeCellType = compositeCellType, probeSelect = probeSelect) {
    splitit <- function(x) {
        split(seq(along=x), x)
    }
    
    p <- getBeta(mSet)
    pd <- as.data.frame(pData(mSet))
    if(is.null(cellTypes)) {
        t <- unique(pd$CellType)
    } else {
        if(!all(cellTypes %in% pd$CellType))
            stop("elements of argument 'cellTypes' is not part of 'mSet$CellType'")
        keep <- which(pd$CellType %in% cellTypes)
        pd <- pd[keep,]
        p <- p[,keep]
    }
    ## make cell type a factor 
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
                                        # get fstats
    ffComp <- rowFtests(p, pd$CellType)
    prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[,i]))
    r <- matrixStats::rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2,-1,0) + ncol(compTable)] <- c("low", "high", "range") 
    tIndexes <- splitit(pd$CellType)
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0,ncol(p))
        x[i] <- 1
        return(rowttests(p, factor(x)))
    })
    
	if (probeSelect == "any"){
		probeList <- lapply(tstatList, function(x) {
			y <- x[x[,"p.value"] < 1e-8,]
			yAny <- y[order(abs(y[,"dm"]), decreasing=TRUE),]      
			c(rownames(yAny)[1:(numProbes*2)])
		})
	} else {
	probeList <- lapply(tstatList, function(x) {
		y <- x[x[,"p.value"] < 1e-8,]
		yUp <- y[order(y[,"dm"], decreasing=TRUE),]
		yDown <- y[order(y[,"dm"], decreasing=FALSE),]
		c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
		})
	}
    
    trainingProbes <- unique(unlist(probeList))
    p <- p[trainingProbes,]
    
    pMeans <- colMeans(p)
    names(pMeans) <- pd$CellType
    
    form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse="+")))
    phenoDF <- as.data.frame(model.matrix(~pd$CellType-1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    if(ncol(phenoDF) == 2) { # two group solution
		X <- as.matrix(phenoDF)
		coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
	} else { # > 2 group solution
		tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
		coefEsts <- tmp$coefEsts
	}
    
    out <- list(coefEsts = coefEsts, compTable = compTable,
                sampleMeans = pMeans)
    return(out)
}

projectCellType <- function(Y, coefCellType, contrastCellType=NULL, nonnegative=TRUE, lessThanOne=FALSE){ 
    if(is.null(contrastCellType))
        Xmat <- coefCellType
    else
        Xmat <- tcrossprod(coefCellType, contrastCellType) 
    
    nCol <- dim(Xmat)[2]
    if(nCol == 2) {
        Dmat <- crossprod(Xmat)
        mixCoef <- t(apply(Y, 2, function(x) { solve(Dmat, crossprod(Xmat, x)) }))
        colnames(mixCoef) <- colnames(Xmat)
        return(mixCoef)
    } else {
        nSubj <- dim(Y)[2]
        
        mixCoef <- matrix(0, nSubj, nCol)
        rownames(mixCoef) <- colnames(Y)
        colnames(mixCoef) <- colnames(Xmat)
        
        if(nonnegative){
            if(lessThanOne) {
                Amat <- cbind(rep(-1, nCol), diag(nCol))
                b0vec <- c(-1, rep(0, nCol))
            } else {
                Amat <- diag(nCol)
                b0vec <- rep(0, nCol)
            }
            for(i in 1:nSubj) {
                obs <- which(!is.na(Y[,i])) 
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve.QP(Dmat, crossprod(Xmat[obs,], Y[obs,i]), Amat, b0vec)$sol
            }
        } else {
            for(i in 1:nSubj) {
                obs <- which(!is.na(Y[,i])) 
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
            }
        }
        return(mixCoef)
    }
}

validationCellType <- function(Y, pheno, modelFix, modelBatch=NULL,
                               L.forFstat = NULL, verbose = FALSE){
    N <- dim(pheno)[1]
    pheno$y <- rep(0, N)
    xTest <- model.matrix(modelFix, pheno)
    sizeModel <- dim(xTest)[2]
    M <- dim(Y)[1]

    if(is.null(L.forFstat)) {
        L.forFstat <- diag(sizeModel)[-1,] # All non-intercept coefficients
        colnames(L.forFstat) <- colnames(xTest) 
        rownames(L.forFstat) <- colnames(xTest)[-1] 
    }

    ## Initialize various containers
    sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
    coefEsts <- matrix(NA, M, sizeModel)
    coefVcovs <- list()

    if(verbose)
        cat("[validationCellType] ")
    for(j in 1:M) { # For each CpG
        ## Remove missing methylation values
        ii <- !is.na(Y[j,])
        nObserved[j] <- sum(ii)
        pheno$y <- Y[j,]

        if(j%%round(M/10)==0 && verbose)
            cat(".") # Report progress

        try({ # Try to fit a mixed model to adjust for plate
            if(!is.null(modelBatch)) {
                fit <- try(lme(modelFix, random=modelBatch, data=pheno[ii,]))
                OLS <- inherits(fit,"try-error") # If LME can't be fit, just use OLS
            } else
                OLS <- TRUE

            if(OLS) {
                fit <- lm(modelFix, data=pheno[ii,])
                fitCoef <- fit$coef
                sigmaResid[j] <- summary(fit)$sigma
                sigmaIcept[j] <- 0
                nClusters[j] <- 0
            } else { 
                fitCoef <- fit$coef$fixed
                sigmaResid[j] <- fit$sigma
                sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
                nClusters[j] <- length(fit$coef$random[[1]])
            }
            coefEsts[j,] <- fitCoef
            coefVcovs[[j]] <- vcov(fit)
            
            useCoef <- L.forFstat %*% fitCoef
            useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
            Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef))/sizeModel
        })
    }
    if(verbose)
        cat(" done\n")
    ## Name the rows so that they can be easily matched to the target data set
    rownames(coefEsts) <- rownames(Y)
    colnames(coefEsts) <- names(fitCoef)
    degFree <- nObserved - nClusters - sizeModel + 1

    ## Get P values corresponding to F statistics
    Pval <- 1-pf(Fstat, sizeModel, degFree)
    
    out <- list(coefEsts=coefEsts, coefVcovs=coefVcovs, modelFix=modelFix, modelBatch=modelBatch,
                sigmaIcept=sigmaIcept, sigmaResid=sigmaResid, L.forFstat=L.forFstat, Pval=Pval,
                orderFstat=order(-Fstat), Fstat=Fstat, nClusters=nClusters, nObserved=nObserved,
                degFree=degFree)
    
    out
}