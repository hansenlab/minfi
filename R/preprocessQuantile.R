preprocessQuantile <- function(object, fixOutliers=TRUE,
                               removeBadSamples=FALSE, badSampleCutoff=10.5,
                               quantileNormalize=TRUE, stratified=TRUE,
                               mergeManifest=FALSE, sex=NULL, verbose=TRUE){
    ## We could use [Genomic]MethylSet if the object has been processed with preprocessRaw()
    if(! (is(object, "RGChannelSet") || is(object, "MethylSet") ||
          is(object, "GenomicMethylSet") ))
        stop("object must be of class 'RGChannelSet' or '[Genomic]MethylSet'")
    if( (is(object, "MethylSet") || is(object, "GenomicMethylSet") ) &&
       preprocessMethod(object)["rg.norm"] != "Raw (no normalization or bg correction)")
        warning("preprocessQuantile has only been tested with 'preprocessRaw'")
    if (!is.null(sex))
        sex <- .checkSex(sex)

    if(verbose) cat("[preprocessQuantile] Mapping to genome.\n")
    object <- mapToGenome(object, mergeManifest = mergeManifest)
    
    if(fixOutliers){
        if(verbose) cat("[preprocessQuantile] Fixing outliers.\n")
        object <- fixMethOutliers(object)
    }
    qc <- getQC(object)
    meds <- (qc$uMed + qc$mMed)/2
    keepIndex <- which(meds > badSampleCutoff)
    if(length(keepIndex) == 0 && removeBadSamples) stop("All samples found to be bad")
    if(length(keepIndex) < ncol(object) && removeBadSamples) {
        if(verbose) cat(sprintf("[preprocessQuantile] Found and removed %s bad samples.\n",
                                ncol(object) - length(keepIndex)))
        object <- object[, keepIndex]
    }

    if (is.null(sex)) {
        object <- addSex(object)
        sex <- pData(object)$predictedSex
    }

    xIndex <- which(seqnames(object) == "chrX")
    yIndex <- which(seqnames(object) == "chrY")
    auIndex <- which(seqnames(object) %in% paste0("chr", 1:22))

    if(quantileNormalize){
        if(verbose) cat("[preprocessQuantile] Quantile normalizing.\n")
        if (!stratified) {
            U <- .qnormNotStratified(getUnmeth(object), auIndex, xIndex, yIndex, sex)
            M <- .qnormNotStratified(getMeth(object), auIndex, xIndex, yIndex, sex)
        } else {
            probeType <- getProbeType(object)
            regionType <- getIslandStatus(object)
            regionType[regionType %in% c("Shelf", "OpenSea")] <- "Far"
            U <- .qnormStratified(getUnmeth(object), auIndex,
                                  xIndex, yIndex, sex, probeType, regionType)
            M <- .qnormStratified(getMeth(object), auIndex,
                                  xIndex, yIndex, sex, probeType, regionType)
        }
    } else {
        U <- getUnmeth(object)
        M <- getMeth(object)
    }
    preprocessMethod <- c(mu.norm="preprocessQuantile", preprocessMethod(object))
    out <- GenomicRatioSet(gr=granges(object), Beta=NULL, M = log2(M/U),
                           CN=log2(U+M), pData=pData(object),
                           annotation=annotation(object),
                           preprocessMethod=preprocessMethod)
    return(out)
}


##quantile normalize but X and Y chromsome by sex
.qnormNotStratified <- function(mat, auIndex, xIndex, yIndex, sex=NULL){
    mat[auIndex,] <- preprocessCore::normalize.quantiles(mat[auIndex,])
    if(!is.null(sex)) {
        sexIndexes <- split(1:ncol(mat),sex)
    } else {
        sexIndexes <- list(U=1:ncol(mat))
    }
    for(i in seq(along=sexIndexes)){
        Index <- sexIndexes[[i]]
        if(length(Index) > 1){
            mat[c(xIndex, yIndex), Index] <- preprocessCore::normalize.quantiles(mat[c(xIndex, yIndex),Index])
        } else{
            warning(sprintf("Only one sample of sex: %s. Not normalizing the sex chromosomes for that sample.",
                            names(sexIndexes)[i]))
        }
    }
    return(mat)
}

.qnormStratified <- function(mat, auIndex, xIndex, yIndex, sex=NULL, probeType, regionType){
    mat[auIndex,] <- .qnormStratifiedHelper(mat[auIndex,], probeType[auIndex], regionType[auIndex])
    if(!is.null(sex)) {
        sexIndexes <- split(1:ncol(mat), sex)
    } else {
        ## If sex if not given, we will assume all samples have same sex
        sexIndexes <- list(1:ncol(mat))
    }
    sexIndexes <- sexIndexes[sapply(sexIndexes, length) > 1]
    for(idxes in sexIndexes) {
        mat[c(xIndex, yIndex), idxes] <- .qnormStratifiedHelper(mat[c(xIndex, yIndex), idxes, drop=FALSE],
                                                                probeType[c(xIndex, yIndex)],
                                                                regionType[c(xIndex, yIndex)])
    }
    return(mat)
}

.qnormStratifiedHelper <- function(mat, probeType, regionType) {
    if(ncol(mat) == 1)
        return(mat)
    if(length(probeType) != length(regionType))
        stop("length of 'probeType' and 'regionType' needs to be the same.")
    if(nrow(mat) != length(probeType))
        stop("'mat' needs to have as many rows as entries in 'probeType'")
    regionTypes <- unique(regionType)
    for(i in seq(along=regionTypes)){
        inRegion <- (regionType == regionTypes[i])
        Index1 <- which(inRegion & probeType == "I")
        Index2 <- which(inRegion & probeType == "II")
        mat[Index2,] <- preprocessCore::normalize.quantiles(mat[Index2,])
        target <- approx(seq(along=Index2), sort(mat[Index2,1]),
                         seq(1,length(Index2), length.out=length(Index1)))$y
        mat[Index1,] <- preprocessCore::normalize.quantiles.use.target(mat[Index1,,drop=FALSE],
                                                                       target)
    }
    return(mat)
}
