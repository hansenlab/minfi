estimateCellCounts <- function (RGset, meanPlot = TRUE, cellType = "Blood",
                                verbose=TRUE, ...) {
    require(quadprog)
    platform <- "450k" # FIXME: get platform for RGset
    referencePkg <- sprintf("FlowSorted.%s.%s", cellType, platform)
    if(!require(referencePkg, character.only = TRUE))
        stop(sprintf("Could not find reference data package for cellType '%s' and platform '%s' (inferred package name is '%s')",
                     cellType, platform, referencePkg))
    data(list = referencePkg)
    referenceRGset <- get(referencePkg)

    if(verbose) cat("[estimateCellCounts] Combining Data with Flow Sorted Data.\n")
    RGsetComb <- combine(RGset, referenceRGset)
    newpd <- data.frame(sampleNames = c(sampleNames(RGset), sampleNames(referenceRGset)),
                        studyIndex = rep(c("user", "reference"),
                        times = c(ncol(RGset), ncol(referenceRGset))),
                        stringsAsFactors=FALSE)
    pData(RGsetComb) <- newpd

    if(verbose) cat("[estimateCellCounts] Normalizing Data Together.\n")
    Mset <- preprocessQuantile(RGsetComb, removeBadSamples = FALSE, ...) 
		
    ## Extracts normalized reference data
    referenceMset <- Mset[, Mset$studyIndex == "reference"]
    pData(referenceMset) <- as(pData(referenceRGset), "DataFrame")
    Mset <- Mset[, Mset$studyIndex == "user"]
    pData(Mset) <- as(pData(RGset), "DataFrame")
    
    if(verbose) cat("[estimateCellCounts] Picking Probes for Composition Estimation.\n")
    compData <- pickCompProbes(referenceMset)
    coefs <- compData$coefEsts
    
    if(verbose) cat("[estimateCellCounts] Estimating Composition.\n")
    counts <- projectWBC(getBeta(Mset)[rownames(coefs), ], coefs)
    rownames(counts) <- sampleNames(RGset)

    if (meanPlot) {
        smeans <- compData$sampleMeans
        smeans <- smeans[order(names(smeans))]
        sampleMeans <- c(colMeans(getBeta(Mset)[rownames(coefs), ]), smeans)

        sampleColors <- c(rep(1, ncol(Mset)), 1 + as.numeric(factor(names(smeans))))
        plot(sampleMeans, pch = 21, bg = sampleColors)
        legend("bottomleft", c("blood", levels(factor(names(smeans)))),
               col = 1:7, pch = 15)
    }
    list(counts = counts, compTable = compData$compTable,
         sortedData = referenceMset)
}

pickCompProbes <- function(Mset, numProbes = 50) {
    splitit <- function(x) {
        split(seq(along=x), x)
    }
    
    p <- getBeta(Mset)
    pd <- as.data.frame(pData(Mset))
    
    ## only keep 6 components from kere
    keep <- which(pd$CellType %in% c("Mono", "Bcell", 
                                     "Gran", "CD4T", "CD8T", "NK"))
    pd <- pd[keep,]
    p <- p[,keep]
    
    ## make cell type a factor 
    pd$CellType <- factor(pd$CellType, 
                          levels = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"))
    
                                        # get fstats
    ffComp <- rowFtests(p, pd$CellType)
    prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[,i]))
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
    names(compTable)[c(1, 9:11)] <- c("Fstat", "low", "high", "range")
    
                                        # t-test by cell type
    tIndexes <- splitit(pd$CellType)
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0,ncol(p))
        x[i] <- 1
        return(rowttests(p, factor(x)))
    })
    
    ## take N up and N down
    probeList <- lapply(tstatList, function(x) {
        y <- x[x[,"p.value"] < 1e-8,]
        yUp <- y[order(y[,"dm"], decreasing=TRUE),]
        yDown <- y[order(y[,"dm"], decreasing=FALSE),]
        c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
    })
    
    trainingProbes <- unlist(probeList)
    p <- p[trainingProbes,]
    
    pMeans <- colMeans(p)
    names(pMeans) <- pd$CellType
    
    mod <- model.matrix(~pd$CellType-1)
    colnames(mod) <- levels(pd$CellType)
    form <- as.formula(sprintf("y ~ %s - 1", colnames(mod), collapse="+"))
    
    tmp <- validationWBC(p,data.frame(mod),form)
    coefEsts <- tmp$coefEsts
    
    out <- list(coefEsts = coefEsts, compTable = compTable,
                sampleMeans = pMeans)
    return(out)
}

projectWBC <- function(Y, coefWBC, contrastWBC=NULL, nonnegative=TRUE, lessThanOne=FALSE){ 
    if(is.null(contrastWBC))
        Xmat <- coefWBC
    else
        Xmat <- coefWBC %*% t(contrastWBC) 

    nCol <- dim(Xmat)[2]
    nSubj <- dim(Y)[2]

    mixCoef <- matrix(0, nSubj, nCol)
    rownames(mixCoef) <- colnames(Y)
    colnames(mixCoef) <- colnames(Xmat)
    
    if(nonnegative){
        library(quadprog)
        
        if(lessThanOne){
            Amat <- cbind(rep(-1,nCol), diag(nCol))
            b0vec <- c(-1,rep(0,nCol))
        } else {
            Amat <- diag(nCol)
            b0vec <- rep(0,nCol)
        }
        
        for(i in 1:nSubj) {
            obs <- which(!is.na(Y[,i])) 
            Dmat <- t(Xmat[obs,]) %*% Xmat[obs,]
            mixCoef[i,] <- solve.QP(Dmat, t(Xmat[obs,]) %*% Y[obs,i], Amat, b0vec)$sol
        }
    } else {
        for(i in 1:nSubj) {
            obs <- which(!is.na(Y[,i])) 
            Dmat <- t(Xmat[obs,]) %*% Xmat[obs,]
            mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
        }
    }
    return(mixCoef)
}

validationWBC <- function(Y, pheno, modelFix, modelBatch=NULL, L.forFstat = NULL){
    N <- dim(pheno)[1]
    pheno$y <- rep(0, N)
    xTest <- model.matrix(modelFix, pheno)
    sizeModel <- dim(xTest)[2]
    M <- dim(Y)[1]

    if(is.null(L.forFstat)) {
        L.forFstat <- diag(sizeModel)[-1,]  #All non-intercept coefficients
        colnames(L.forFstat) <- colnames(xTest) 
        rownames(L.forFstat) <- colnames(xTest)[-1] 
    }

    ## Initialize various containers
    sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
    coefEsts <- matrix(NA, M, sizeModel)
    coefVcovs <- list()

    for(j in 1:M) { # For each CpG
        ## Remove missing methylation values
        ii <- !is.na(Y[j,])
        nObserved[j] <- sum(ii)
        pheno$y <- Y[j,]

        if(j%%round(M/10)==0)
            cat(j,"\n") # Report progress

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

