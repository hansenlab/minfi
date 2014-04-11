getSubset <- function(counts, subset){
    x <- numeric(0)
    for(i in 1:3){
        x <- c(x,sample(seq(1, length(counts), by = 1)[counts == i], subset))
    }
    return(seq(1, length(counts)) %in% x)
}

bgIntensitySwan <- function(rgSet){
    grnMed <- colMedians(getGreen(rgSet)[getControlAddress(rgSet, controlType = "NEGATIVE"), ])
    redMed <- colMedians(getRed(rgSet)[getControlAddress(rgSet, controlType = "NEGATIVE"), ])
    return(rowMeans(cbind(grnMed, redMed)))
}

preprocessSWAN <- function(rgSet, mSet = NULL, verbose = FALSE){
    if(is.null(mSet))
        mSet <- preprocessRaw(rgSet)
    typeI <- getProbeInfo(rgSet, type = "I")[, c("Name", "nCpG")]
    typeII <- getProbeInfo(rgSet, type = "II")[, c("Name", "nCpG")]
    ## This next part should be fixed so it becomes more elegant
    CpG.counts <- rbind(typeI, typeII)
    CpG.counts$Name <- as.character(CpG.counts$Name)
    CpG.counts$Type <- rep(c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
    names(CpG.counts)[2] <- "CpGs"
    counts <- CpG.counts[CpG.counts$Name %in% featureNames(mSet),]
    subset <- min(table(counts$CpGs[counts$Type == "I" & counts$CpGs %in% 1:3]),
                  table(counts$CpGs[counts$Type == "II" & counts$CpGs %in% 1:3]))
    bg <- bgIntensitySwan(rgSet)
    methData <- getMeth(mSet)
    unmethData <- getUnmeth(mSet)
    xNormSet <- vector("list", 2)
    xNormSet[[1]] <- getSubset(counts$CpGs[counts$Type=="I"], subset)
    xNormSet[[2]] <- getSubset(counts$CpGs[counts$Type=="II"], subset)
    normMethData <- matrix(NA_real_, ncol = ncol(methData), nrow = nrow(methData))
    colnames(normMethData) <- sampleNames(mSet)
    normUnmethData <- normMethData
    normSet <- mSet
    for(i in 1:ncol(mSet)) {
        if(verbose) cat(sprintf("[preprocessSwan] Normalizing array %d of %d\n", i, ncol(mSet)))
        normMeth <- normaliseChannel(methData[rownames(methData) %in% counts$Name[counts$Type=="I"], i],
                                     methData[rownames(methData) %in% counts$Name[counts$Type=="II"], i],
                                     xNormSet, bg[i])
        normMethData[, i] <- normMeth
        normUnmeth <- normaliseChannel(unmethData[rownames(unmethData) %in% counts$Name[counts$Type=="I"], i],
                                       unmethData[rownames(unmethData) %in% counts$Name[counts$Type=="II"], i],
                                       xNormSet, bg[i])
        normUnmethData[, i] <- normUnmeth
    }
    rownames(normMethData) <- names(normMeth)
    rownames(normUnmethData) <- names(normUnmeth)
    assayDataElement(normSet, "Meth") <- normMethData
    assayDataElement(normSet, "Unmeth") <- normUnmethData
    normSet@preprocessMethod <- c(rg.norm = sprintf("SWAN (based on a MethylSet preprocesses as '%s'",
                                          preprocessMethod(mSet)[1]),
                                  minfi = as.character(packageVersion("minfi")),
                                  manifest = as.character(packageVersion("IlluminaHumanMethylation450kmanifest")))
    normSet
}

normaliseChannel <- function(intensityI, intensityII, xNormSet, bg) {
    xTarget <- aveQuantile(list(intensityI[xNormSet[[1]]], intensityII[xNormSet[[2]]]))
    xNorm <- unlist(subsetQuantileNorm(list(intensityI, intensityII), xNormSet, xTarget, bg))
    names(xNorm) <- c(names(intensityI), names(intensityII))
    xNorm
}

aveQuantile <- function(X) {
    nbrOfChannels <- length(X)
    if (nbrOfChannels == 1) {
        return(X)
    }
    nbrOfObservations <- unlist(lapply(X, FUN = length), use.names = FALSE)
    maxNbrOfObservations <- max(nbrOfObservations)
    if (maxNbrOfObservations == 1) {
        return(X)
    }
    ## nbrOfFiniteObservations <- rep(maxNbrOfObservations, times = nbrOfChannels)
    quantiles <- (0:(maxNbrOfObservations - 1))/(maxNbrOfObservations - 1)
    xTarget <- vector("double", maxNbrOfObservations)
    for (cc in 1:nbrOfChannels) {
        Xcc <- X[[cc]]
        Scc <- sort(Xcc)
        nobs <- length(Scc)
        if (nobs < maxNbrOfObservations) {
            ## tt <- !is.na(Xcc)
            bins <- (0:(nobs - 1))/(nobs - 1)
            Scc <- approx(x = bins, y = Scc, xout = quantiles,ties = "ordered")$y
        }
        xTarget <- xTarget + Scc
    }
    rm(Scc, Xcc)
    xTarget <- xTarget/nbrOfChannels
    xTarget
}

subsetQuantileNorm <- function(x, xNormSet, xTarget, bg) {
    for(i in 1:length(x)){
        n <- length(x[[i]])
        nTarget <- length(xTarget)
        nNormSet <- sum(xNormSet[[i]])
        
        if(nNormSet != nTarget){
            targetQuantiles <- (0:(nTarget - 1))/(nTarget - 1)
            r <- rank(x[xNormSet[,i], i])
            xNew <-(r - 1)/(nNormSet - 1)
            xNew <- xNew[order(xNew)]
            xNorm <- approx(x = targetQuantiles, y = xTarget, xout = xNew, ties = "ordered", rule = 2)$y
        } else {
            xNorm<-xTarget
        }
        
        r <- rank(x[[i]])
        xNew <-(r - 1)/(n - 1)
        quantiles <- xNew[xNormSet[[i]]]
        quantiles <- quantiles[order(quantiles)]
        xmin <- min(x[[i]][xNormSet[[i]]]) #get min value from subset
        xmax <- max(x[[i]][xNormSet[[i]]]) #get max value from subset
        kmax <- which(xNew > max(quantiles)) 
        kmin<- which(xNew < min(quantiles))
        offsets.max <- x[[i]][kmax]-xmax
        offsets.min <- x[[i]][kmin]-xmin
        x[[i]] <- approx(x = quantiles, y = xNorm, xout = xNew, ties = "ordered")$y #interpolate
        x[[i]][kmax]<- max(xNorm) + offsets.max
        x[[i]][kmin]<- min(xNorm) + offsets.min
        x[[i]] = ifelse(x[[i]] <= 0, bg, x[[i]])    
    }  
    x
}


