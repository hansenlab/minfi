getSubset <- function(counts, subset){
    x <- integer(0)
    for (i in 1:3) {
        x <- c(x, sample(seq.int(1, length(counts))[counts == i], subset))
    }
    return(seq.int(1, length(counts)) %in% x)
}

bgIntensitySwan <- function(rgSet) {
    rows <- match(getControlAddress(rgSet, controlType = "NEGATIVE"),
                  rownames(rgSet))
    grnMed <- colMedians(
        x = getGreen(rgSet),
        rows = rows)
    redMed <- colMedians(
        x = getRed(rgSet),
        rows = rows)
    return(rowMeans(cbind(grnMed, redMed)))
}

# x is either `Meth` or `Unmeth`
setGeneric(
    ".preprocessSWAN",
    function(x, xNormSet, counts, bg) standardGeneric(".preprocessSWAN"),
    signature = "x")

# TODO: Parallelise with BiocParallel
setMethod(".preprocessSWAN", "matrix", function(x, xNormSet, counts, bg) {
    # NOTE: SWAN can return non-integer values, so fill a numeric matrix
    normalized_x <- matrix(NA_real_,
                           ncol = ncol(x),
                           nrow = nrow(x),
                           dimnames = dimnames(x))
    typeI_idx <- rownames(x) %in% counts$Name[counts$Type == "I"]
    typeII_idx <- rownames(x) %in% counts$Name[counts$Type == "II"]
    for (j in seq_len(ncol(x))) {
        normalized_x[, j] <- normaliseChannel(
            intensityI = x[typeI_idx, j],
            intensityII = x[typeII_idx, j],
            xNormSet = xNormSet,
            bg = bg[j])
    }
    normalized_x
})

# TODO: Parallelise with BiocParallel
setMethod(".preprocessSWAN", "DelayedMatrix", function(x, xNormSet, counts, bg) {
    # Set up intermediate RealizationSink objects of appropriate dimensions and
    # type
    # NOTE: This is ultimately coerced to the output DelayedMatrix object,
    #       `normalized_x`
    # NOTE: SWAN can return non-integer values, so fill a numeric sink
    normalized_x_type <- "double"
    normalized_x_sink <- DelayedArray:::RealizationSink(
        dim = c(nrow(x), ncol(x)),
        dimnames = dimnames(x),
        type = normalized_x_type)

    # Set up ArrayGrid instances over `x`.
    # We're going to walk over columns of `x` and write to columns
    # of `normalized_x_sink`. This imposes some requirement(s):
    # 1. Need to increase the block length so each block is made of at least
    #    one column of `x`
    max_block_len <- max(
        nrow(x),
        DelayedArray:::get_max_block_length(normalized_x_type))
    x_grid <- defaultGrid(x, max_block_len)

    blockApplyWithRealization(
        x = x,
        xNormSet = xNormSet,
        counts = counts,
        bg = bg,
        FUN = .preprocessSWAN,
        grid = x_grid,
        sink = normalized_x_sink,
        BPREDO = list(),
        BPPARAM = SerialParam())

    as(normalized_x_sink, "DelayedArray")
})

# TODO: Parallelise with BiocParallel
preprocessSWAN <- function(rgSet, mSet = NULL, verbose = FALSE){
    if (is.null(mSet)) {
        MSet <- preprocessRaw(rgSet)
    }
    typeI <- getProbeInfo(rgSet, type = "I")[, c("Name", "nCpG")]
    typeII <- getProbeInfo(rgSet, type = "II")[, c("Name", "nCpG")]
    ## This next part should be fixed so it becomes more elegant
    CpG.counts <- rbind(typeI, typeII)
    # TODO: Unclear why this is necessary
    CpG.counts$Name <- as.character(CpG.counts$Name)
    CpG.counts$Type <- rep(c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
    counts <- CpG.counts[CpG.counts$Name %in% rownames(MSet), ]
    subset <- min(
        table(counts$nCpG[counts$Type == "I" & counts$nCpG %in% 1:3]),
        table(counts$nCpG[counts$Type == "II" & counts$nCpG %in% 1:3]))
    bg <- bgIntensitySwan(rgSet)
    Meth <- getMeth(MSet)
    Unmeth <- getUnmeth(MSet)
    xNormSet <- lapply(c("I", "II"), function(type) {
        getSubset(counts$nCpG[counts$Type == type], subset)
    })
    assay(MSet, "Meth") <- .preprocessSWAN(
        x = Meth,
        xNormSet = xNormSet,
        counts = counts,
        bg = bg)
    assay(MSet, "Unmeth") <- .preprocessSWAN(
        x = Unmeth,
        xNormSet = xNormSet,
        counts = counts,
        bg = bg)
    MSet@preprocessMethod <- c(
        rg.norm = sprintf("SWAN (based on a MethylSet preprocesses as '%s'",
                          preprocessMethod(MSet)[1]),
        minfi = as.character(packageVersion("minfi")),
        manifest = as.character(
            packageVersion(.getManifestString(annotation(rgSet)))))
    MSet
}

# TODO: Profile
normaliseChannel <- function(intensityI, intensityII, xNormSet, bg) {
    xTarget <- aveQuantile(
        list(intensityI[xNormSet[[1]]], intensityII[xNormSet[[2]]]))
    xNorm <- unlist(
        subsetQuantileNorm(
            list(intensityI, intensityII), xNormSet, xTarget, bg))
    names(xNorm) <- c(names(intensityI), names(intensityII))
    xNorm
}

# TODO: Profile
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

# TODO: Profile
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

# ------------------------------------------------------------------------------
# TODOs
#

# TODO: Choose more informative names for variables
