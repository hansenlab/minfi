
.detectionP_matrix <- function(dimnames, Red, Green, controlIdx, TypeI.Red, TypeI.Green, TypeII) {
    # Compute summary statistics needed for calculations
    rBg <- Red[controlIdx, , drop = FALSE]
    rMu <- colMedians(rBg)
    rSd <- colMads(rBg)
    gBg <- Green[controlIdx, , drop = FALSE]
    gMu <- colMedians(gBg)
    gSd <- colMads(gBg)

    # Set up output matrix with appropriate dimensions and type
    dim <- lengths(dimnames)
    detP <- matrix(NA_real_,
                   nrow = dim[[1L]],
                   ncol = dim[[2L]],
                   dimnames = dimnames)

    # Fill output matrix
    for (j in seq_len(ncol(detP))) {
        # Type I Red
        intensity <- Red[TypeI.Red$AddressA, j] + Red[TypeI.Red$AddressB, j]
        detP[TypeI.Red$Name, j] <- 1 - pnorm(intensity, mean = rMu[j] * 2, sd = rSd[j] * 2)
        # Type I Green
        intensity <- Green[TypeI.Green$AddressA, j] + Green[TypeI.Green$AddressB, j]
        detP[TypeI.Green$Name, j] <- 1 - pnorm(intensity, mean = gMu[j] * 2, sd = gSd[j] * 2)
        # Type II
        intensity <- Red[TypeII$AddressA, j] + Green[TypeII$AddressA, j]
        detP[TypeII$Name, j] <- 1 - pnorm(intensity, mean = rMu[j] + gMu[j], sd = rSd[j] + gSd[j])
    }

    # Return output matrix
    detP
}

.detectionP_DelayedMatrix <- function(dimnames, Red, Green, controlIdx, TypeI.Red, TypeI.Green,
                                      TypeII, BACKEND = getRealizationBackend()) {
    dim <- lengths(dimnames)
    # TODO: 'sink' doesn't have conformable dimensions to 'Red' or 'Green', so
    #       DelayedArray:::block_MAPPLY() will refuse to write to it
    sink <- DelayedArray:::RealizationSink(dim = dim,
                                           dimnames = dimnames,
                                           type = "double")
    # We're going to walk over columns of 'Red' and 'Green', which mean we need to increase
    # the block length so each block is made of at least one column.
    max_block_len <- max(DelayedArray:::get_max_block_length(DelayedArray::type(Red)),
                         dim[[1L]])
    detP_columns <- DelayedArray:::block_MAPPLY(function(Red, Green) {
        # Compute detP on block
        .detectionP_matrix(dimnames = list(dimnames[[1L]], colnames(Red)),
                           Red = as.matrix(Red),
                           Green = as.matrix(Green),
                           controlIdx = controlIdx,
                           TypeI.Red = TypeI.Red,
                           TypeI.Green = TypeI.Green,
                           TypeII = TypeII)
    }, Red = Red, Green = Green, sink = sink, max_block_len = max_block_len)

    # Return output matrix
    do.call(cbind, detP_columns)
}

# TODO: In the release version of minfi, detectionP() always returns a matrix. This
#       version can return a matrix or a DelayedMatrix. Is this a good idea and/or
#       necessary?
detectionP <- function(rgSet, type = "m+u") {
    .isRGOrStop(rgSet)

    # Extract data to pass to low-level function that constructs 'detP'
    locusNames <- getManifestInfo(rgSet, "locusNames")
    sampleNames <- sampleNames(rgSet)
    dimnames <- list(locusNames, sampleNames)
    controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")
    Red <- getRed(rgSet)
    Green <- getGreen(rgSet)
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    TypeII <- getProbeInfo(rgSet, type = "II")

    # Construct 'detP'
    if (is(Red, "DelayedMatrix") || is(Green, "DelayedMatrix")) {
        detP <- .detectionP_DelayedMatrix(dimnames, Red, Green, TypeI.Red, TypeI.Green,
                                          TypeII)
    } else if (is.matrix(Red) && is.matrix(Green)) {
        detP <- .detectionP_matrix()
    }
    detP
}
