# Internal functions -----------------------------------------------------------

getSubset <- function(counts, subset){
    x <- integer(0)
    for (i in 1:3) {
        x <- c(x, sample(seq.int(1, length(counts))[counts == i], subset))
    }
    seq.int(1, length(counts)) %in% x
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
    rowMeans2(cbind(grnMed, redMed))
}

# TODO: Profile and tidy
normaliseChannel <- function(intensityI, intensityII, xNormSet, bg) {
    xTarget <- aveQuantile(
        list(intensityI[xNormSet[[1]]], intensityII[xNormSet[[2]]]))
    xNorm <- unlist(
        subsetQuantileNorm(
            list(intensityI, intensityII), xNormSet, xTarget, bg))
    names(xNorm) <- c(names(intensityI), names(intensityII))
    xNorm
}

# TODO: Profile and tidy
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

# TODO: Profile and tidy
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

# Internal generics ------------------------------------------------------------

# `x` is either `Meth` or `Unmeth`
# `...` are additional arguments passed to methods.
setGeneric(
    ".preprocessSWAN",
    function(x, xNormSet, counts, bg, ...) standardGeneric(".preprocessSWAN"),
    signature = "x")

# Internal methods -------------------------------------------------------------

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

setMethod(
    ".preprocessSWAN",
    "DelayedMatrix",
    function(x, xNormSet, counts, bg, BPREDO = list(),
             BPPARAM = SerialParam()) {
        # Set up intermediate RealizationSink object of appropriate dimensions
        # and type
        # NOTE: This is ultimately coerced to the output DelayedMatrix object
        # NOTE: SWAN can return non-integer values, so fill a "double" sink
        ans_type <- "double"
        sink <- DelayedArray::AutoRealizationSink(
            dim = c(nrow(x), ncol(x)),
            dimnames = dimnames(x),
            type = ans_type)
        on.exit(close(sink))

        # Coerce `bg` to a row vector
        # NOTE: This is required in order to iterate in "parallel" over `x`,
        #       `bg`, and `sink`.
        bg <- matrix(bg, ncol = length(bg))

        # Set up ArrayGrid instances over `x` as well as "parallel" ArrayGrid
        # instances over `bg` and `sink`.
        x_grid <- colAutoGrid(x)
        bg_grid <- RegularArrayGrid(
            refdim = dim(bg),
            spacings = c(nrow(bg), ncol(bg) / length(x_grid)))
        sink_grid <- RegularArrayGrid(
            refdim = dim(sink),
            spacings = c(nrow(sink), ncol(sink) / length(x_grid)))
        # Sanity check ArrayGrid objects have the same dim
        stopifnot(dim(x_grid) == dim(sink_grid))

        # Loop over column-blocks of `x`, perform SWAN normalization, and write
        # to `normalized_x_sink`
        blockMapplyWithRealization(
            FUN = .preprocessSWAN,
            x = x,
            bg = bg,
            MoreArgs = list(
                xNormSet = xNormSet,
                counts = counts),
            sinks = list(sink),
            dots_grids = list(x_grid, bg_grid),
            sinks_grids = list(sink_grid),
            BPREDO = BPREDO,
            BPPARAM = BPPARAM)

        # Return as DelayedMatrix
        as(sink, "DelayedArray")
    })

# Exported functions -----------------------------------------------------------

preprocessSWAN <- function(rgSet, mSet = NULL, verbose = FALSE) {
    .isRGOrStop(rgSet)

    # Extract data to pass to low-level functions that construct normalized `M`
    # and `U`
    if (is.null(mSet)) {
        mSet <- preprocessRaw(rgSet)
    } else {
        .isMethylOrStop(mSet)
    }
    typeI <- getProbeInfo(rgSet, type = "I")[, c("Name", "nCpG")]
    typeII <- getProbeInfo(rgSet, type = "II")[, c("Name", "nCpG")]
    # TODO This next part should be fixed so it becomes more elegant
    CpG.counts <- rbind(typeI, typeII)
    # TODO: Unclear why this is necessary
    CpG.counts$Name <- as.character(CpG.counts$Name)
    CpG.counts$Type <- rep(c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
    counts <- CpG.counts[CpG.counts$Name %in% rownames(mSet), ]
    subset <- min(
        table(counts$nCpG[counts$Type == "I" & counts$nCpG %in% 1:3]),
        table(counts$nCpG[counts$Type == "II" & counts$nCpG %in% 1:3]))
    bg <- bgIntensitySwan(rgSet)
    Meth <- getMeth(mSet)
    Unmeth <- getUnmeth(mSet)
    xNormSet <- lapply(c("I", "II"), function(type) {
        getSubset(counts$nCpG[counts$Type == type], subset)
    })

    # Construct normalized data
    M <- .preprocessSWAN(
        x = Meth,
        xNormSet = xNormSet,
        counts = counts,
        bg = bg)
    U <- .preprocessSWAN(
        x = Unmeth,
        xNormSet = xNormSet,
        counts = counts,
        bg = bg)

    # Construct MethylSet
    assay(mSet, "Meth") <- M
    assay(mSet, "Unmeth") <- U
    mSet@preprocessMethod <- c(
        rg.norm = sprintf("SWAN (based on a MethylSet preprocesses as '%s'",
                          preprocessMethod(mSet)[1]),
        minfi = as.character(packageVersion("minfi")),
        manifest = as.character(
            packageVersion(.getManifestString(annotation(rgSet)))))
    mSet
}

# TODOs ------------------------------------------------------------------------

# TODO: Choose more informative names for variables
