# Internal functions -----------------------------------------------------------

.qnormStratifiedHelper <- function(mat, probeType, regionType) {
    # Check inputs
    if (ncol(mat) == 1) return(mat)
    if (length(probeType) != length(regionType)) {
        stop("length of 'probeType' and 'regionType' needs to be the same.")
    }
    if (nrow(mat) != length(probeType)) {
        stop("'mat' needs to have as many rows as entries in 'probeType'")
    }

    # Quantile normalize probes in each region type
    regionTypes <- unique(regionType)
    for (i in seq_along(regionTypes)) {
        inRegion <- (regionType == regionTypes[i])
        Index1 <- which(inRegion & probeType == "I")
        Index2 <- which(inRegion & probeType == "II")
        mat[Index2,] <- normalize.quantiles(mat[Index2, ])
        # NOTE: The following code is easy to understand, but we are not using
        #       it because we want the results to stay stable. It gives almost
        #       the same results as the code below:
        # mat[Index1,] <- normalize.quantiles.use.target(
        #     x = mat[Index1, , drop = FALSE],
        #     target = mat[Index2, 1, drop = TRUE])
        target <- approx(
            x = seq(along = Index2),
            y = sort(mat[Index2, 1]),
            xout = seq(1, length(Index2), length.out = length(Index1)))$y
        mat[Index1, ] <- normalize.quantiles.use.target(
            x = mat[Index1, , drop = FALSE],
            target = target)
    }

    mat
}

.qnormStratified <- function(mat, auIndex, xIndex, yIndex, sex = NULL,
                             probeType, regionType) {
    # Normalize probes on the autosomes
    mat[auIndex, ] <- .qnormStratifiedHelper(
        mat = mat[auIndex, ],
        probeType = probeType[auIndex],
        regionType = regionType[auIndex])

    # Normalize probes on the sex chromosomes
    if (is.null(sex)) {
        # NOTE: If sex if not given, we will assume all samples have same sex
        sexIndexes <- list(seq_len(ncol(mat)))
    } else {
        sexIndexes <- split(seq_len(ncol(mat)), sex)

    }
    sexIndexes <- sexIndexes[lengths(sexIndexes) > 1]
    for (idxes in sexIndexes) {
        mat[c(xIndex, yIndex), idxes] <- .qnormStratifiedHelper(
            mat = mat[c(xIndex, yIndex), idxes, drop = FALSE],
            probeType = probeType[c(xIndex, yIndex)],
            regionType = regionType[c(xIndex, yIndex)])
    }

    mat
}

# Quantile normalize but do chrX and chrY separately by sex
.qnormNotStratified <- function(mat, auIndex, xIndex, yIndex, sex = NULL) {
    # Normalize probes on the autosomes
    mat[auIndex,] <- normalize.quantiles(mat[auIndex, ])

    # Normalize probes on the sex chromosomes
    if (is.null(sex)) {
        sexIndexes <- list(U = seq_len(ncol(mat)))
    } else {
        sexIndexes <- split(seq_len(ncol(mat)), sex)
    }
    for (i in seq_along(sexIndexes)) {
        Index <- sexIndexes[[i]]
        if (length(Index) > 1) {
            mat[c(xIndex, yIndex), Index] <- normalize.quantiles(
                x = mat[c(xIndex, yIndex), Index])
        } else {
            warning(
                sprintf("Only one sample of sex: %s. Not normalizing the sex chromosomes for that sample.",
                        names(sexIndexes)[i]))
        }
    }

    mat
}

# Exported functions -----------------------------------------------------------

preprocessQuantile <- function(object, fixOutliers = TRUE,
                               removeBadSamples = FALSE, badSampleCutoff = 10.5,
                               quantileNormalize = TRUE, stratified = TRUE,
                               mergeManifest = FALSE, sex = NULL,
                               verbose = TRUE) {
    # Check inputs
    .isMatrixBackedOrStop(object, "preprocessQuantile")
    # NOTE (Kasper): We could use [Genomic]MethylSet if the object has been
    #                processed with preprocessRaw()
    # TODO: Add the above support?
    if (!(is(object, "RGChannelSet") ||
          is(object, "MethylSet") ||
          is(object, "GenomicMethylSet"))) {
        stop("object must be of class 'RGChannelSet' or '[Genomic]MethylSet'")
    }
    if ((is(object, "MethylSet") || is(object, "GenomicMethylSet")) &&
        (is.na(preprocessMethod(object)["rg.norm"]) ||
         preprocessMethod(object)["rg.norm"] !=
         "Raw (no normalization or bg correction)")) {
        warning("preprocessQuantile has only been tested with 'preprocessRaw'")
    }
    if (.is27k(object) && stratified) {
        stratified <- FALSE
        warning("The stratification option is not available for 27k arrays.")
    }

    # Map to genome
    if (verbose) message("[preprocessQuantile] Mapping to genome.")
    object <- mapToGenome(object, mergeManifest = mergeManifest)

    # Get sex
    if (is.null(sex)) {
        object <- addSex(object)
        sex <- colData(object)$predictedSex
    } else {
        sex <- .checkSex(sex)
    }

    # Fix outliers
    if (fixOutliers) {
        if (verbose) message("[preprocessQuantile] Fixing outliers.")
        object <- fixMethOutliers(object)
    }

    # Run QC
    qc <- getQC(object)
    meds <- (qc$uMed + qc$mMed) / 2
    keepIndex <- which(meds > badSampleCutoff)
    if (length(keepIndex) == 0 && removeBadSamples) {
        stop("All samples found to be bad")
    }
    if (length(keepIndex) < ncol(object) && removeBadSamples) {
        if (verbose) {
            message(
                sprintf("[preprocessQuantile] Found and removed %s bad samples",
                        ncol(object) - length(keepIndex)))
        }
        object <- object[, keepIndex]
    }

    xIndex <- which(seqnames(object) == "chrX")
    yIndex <- which(seqnames(object) == "chrY")
    auIndex <- which(seqnames(object) %in% paste0("chr", 1:22))

    # Quantile normalize Meth and Unmeth
    U <- getUnmeth(object)
    M <- getMeth(object)
    if (quantileNormalize) {
        if (verbose) message("[preprocessQuantile] Quantile normalizing.")
        if (!stratified) {
            U <- .qnormNotStratified(
                mat = U,
                auIndex = auIndex,
                xIndex = xIndex,
                yIndex = yIndex,
                sex = sex)
            if (!is.null(getRealizationBackend())) {
                U <- realize(U)
            }
            M <- .qnormNotStratified(
                mat = M,
                auIndex = auIndex,
                xIndex = xIndex,
                yIndex = yIndex,
                sex = sex)
            if (!is.null(getRealizationBackend())) {
                M <- realize(M)
            }
        } else {
            probeType <- getProbeType(object)
            regionType <- getIslandStatus(object)
            regionType[regionType %in% c("Shelf", "OpenSea")] <- "Far"
            U <- .qnormStratified(
                mat = U,
                auIndex = auIndex,
                xIndex = xIndex,
                yIndex = yIndex,
                sex = sex,
                probeType = probeType,
                regionType = regionType)
            if (!is.null(getRealizationBackend())) {
                U <- realize(U)
            }
            M <- .qnormStratified(
                mat = M,
                auIndex = auIndex,
                xIndex = xIndex,
                yIndex = yIndex,
                sex = sex,
                probeType = probeType,
                regionType = regionType)
            if (!is.null(getRealizationBackend())) {
                M <- realize(M)
            }
        }
    }

    # Construct output GenomicRatioSet
    preprocessMethod <- c(
        mu.norm = "preprocessQuantile",
        preprocessMethod(object))
    GenomicRatioSet(
        gr = granges(object),
        Beta = NULL,
        M = log2(M / U),
        CN = log2(U + M),
        colData = colData(object),
        annotation = annotation(object),
        preprocessMethod = preprocessMethod)
}
