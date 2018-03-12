.preProcessRaw_matrix <- function(type, dimnames, Red, Green, TypeI.Red, TypeI.Green,
                                  TypeII) {

    # Set up output matrices with appropriate dimensions and type
    dim <- lengths(dimnames)
    M <- matrix(.NA_type(type),
                nrow = dim[[1L]],
                ncol = dim[[2L]],
                dimnames = dimnames)
    U <- M

    # Fill output matrices
    M[TypeI.Red$Name, ] <- Red[TypeI.Red$AddressB, ]
    M[TypeI.Green$Name, ] <- Green[TypeI.Green$AddressB, ]
    M[TypeII$Name, ] <- Green[TypeII$AddressA, ]

    U[TypeI.Red$Name, ] <- Red[TypeI.Red$AddressA, ]
    U[TypeI.Green$Name, ] <- Green[TypeI.Green$AddressA, ]
    U[TypeII$Name, ] <- Red[TypeII$AddressA, ]

    # Return output matrices
    list(M = M, U = U)
}

.preProcessRaw_DelayedMatrix <- function(type, dimnames, Red, Green, TypeI.Red,
                                         TypeI.Green, TypeII) {

    # Set up intermediate RealizationSink objects of appropriate dimensions and type
    dim <- lengths(dimnames)
    M_sink <- DelayedArray:::RealizationSink(dim = dim,
                                             dimnames = dimnames,
                                             type = type)
    # NOTE: Don't do `U_sink <- M_sink` or else these will reference the same
    #       object!
    U_sink <- DelayedArray:::RealizationSink(dim = dim,
                                             dimnames = dimnames,
                                             type = type)

    # Fill intermediate RealizationSink objects
    # M
    M_i <- c(TypeI.Red$Name, TypeI.Green$Name, TypeII$Name)
    M_value <- rbind(Red[TypeI.Red$AddressB, ],
                     Green[TypeI.Green$AddressB, ],
                     Green[TypeII$AddressA, ])
    subassignRowsToRealizationSink(M_sink, M_i, M_value)
    # U
    U_i <- c(TypeI.Red$Name, TypeI.Green$Name, TypeII$Name)
    U_value <- rbind(Red[TypeI.Red$AddressA, ],
                     Green[TypeI.Green$AddressA, ],
                     Red[TypeII$AddressA, ])
    subassignRowsToRealizationSink(U_sink, U_i, U_value)

    # Coerce intermediate RealizationSink objects to output DelayedMatrix objects
    M <- as(M_sink, "DelayedArray")
    U <- as(U_sink, "DelayedArray")

    # Return output DelayedMatrix objects
    list(M = M, U = U)
}

# TODO: Add 'BACKEND' argument so that output can have different backend to input?
#       Or should that be driven by setRealizationBackend()?
# TODO: Add BPPARAM (and other BiocParallel args, e.g., BPREDO?)?
preprocessRaw <- function(rgSet) {
    .isRGOrStop(rgSet)

    # Extract data to pass to low-level function that constructs 'M' and 'U'
    locusNames <- getManifestInfo(rgSet, "locusNames")
    sampleNames <- sampleNames(rgSet)
    dimnames <- list(locusNames, sampleNames)
    Red <- getRed(rgSet)
    Green <- getGreen(rgSet)
    type <- .highestType(Red, Green)
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    TypeII <- getProbeInfo(rgSet, type = "II")

    # Construct 'M' and 'U'
    if (is(Red, "DelayedMatrix") || is(Green, "DelayedMatrix")) {
        val <- .preProcessRaw_DelayedMatrix(type, dimnames, Red, Green, TypeI.Red,
                                            TypeI.Green, TypeII)
    } else if (is.matrix(Red) && is.matrix(Green)) {
        val <- .preProcessRaw_matrix(type, dimnames, Red, Green, TypeI.Red,
                                     TypeI.Green, TypeII)
    }

    # Construct MethylSet
    out <- MethylSet(Meth = val$M, Unmeth = val$U, colData = colData(rgSet),
                     annotation = annotation(rgSet), metadata = metadata(rgSet))
    ## TODO:
    ## The manifest package version is currently not updated since 'packageVersion(getManifest(rgSet))' fails.
    ## packageVersion expects a string
    out@preprocessMethod <- c(
        rg.norm = "Raw (no normalization or bg correction)",
        minfi = as.character(packageVersion("minfi")),
        manifest = as.character(packageVersion(.getManifestString(rgSet@annotation))))
    #packageVersion(getManifest(rgSet)))
    out
}

normalize.illumina.control <- function(rgSet, reference=1) {
    ## This function returns an rgset, not a methylset
    ## code duplication
    Green <- getGreen(rgSet)
    Red <- getRed(rgSet)

    if (.is450k(rgSet) || .isEPIC(rgSet)) {
        AT.controls <- getControlAddress(rgSet, controlType = c("NORM_A", "NORM_T"))
        CG.controls <- getControlAddress(rgSet, controlType = c("NORM_C", "NORM_G"))
    }
    if (.is27k(rgSet)) {
        AT.controls <- getControlAddress(rgSet, controlType = "Normalization-Red")
        CG.controls <- getControlAddress(rgSet, controlType = "Normalization-Green")
    }
    Green.avg <- colMeans(Green[CG.controls, , drop = FALSE])
    Red.avg <- colMeans(Red[AT.controls, , drop = FALSE])
    ref <- (Green.avg + Red.avg)[reference]/2
    if(is.na(ref))
        stop("perhaps 'reference' refer to an array that is not present.")
    Green.factor <- ref/Green.avg
    Red.factor <- ref/Red.avg
    # TODO: Need a sweep,DelayedMatrix-method
    #       https://github.com/Bioconductor/DelayedArray/issues/8
    Green <- sweep(Green, 2, FUN = "*", Green.factor)
    Red <- sweep(Red, 2, FUN = "*", Red.factor)
    assay(rgSet, "Green") <- Green
    assay(rgSet, "Red") <- Red
    rgSet
}

bgcorrect.illumina <- function(rgSet) {
    .isRGOrStop(rgSet)
    Green <- getGreen(rgSet)
    Red <- getRed(rgSet)
    if (.is450k(rgSet) || .isEPIC(rgSet)) {
        NegControls <- getControlAddress(rgSet, controlType = "NEGATIVE")
    }
    if (.is27k(rgSet)) {
        NegControls <- getControlAddress(rgSet, controlType = "Negative")
    }
    Green.bg <- apply(Green[NegControls, , drop = FALSE], 2, function(xx) {
        sort(as.vector(xx))[31]
    })
    Red.bg <- apply(Red[NegControls, , drop = FALSE], 2, function(xx) {
        sort(as.vector(xx))[31]
    })
    # TODO: Need a sweep,DelayedMatrix-method
    #       https://github.com/Bioconductor/DelayedArray/issues/8
    Green <- pmax(sweep(Green, 2, Green.bg), 0)
    Red <- pmax(sweep(Red, 2, Red.bg), 0)
    assay(rgSet, "Green") <- Green
    assay(rgSet, "Red") <- Red
    rgSet
}


preprocessIllumina <- function(rgSet, bg.correct = TRUE, normalize = c("controls", "no"),
                               reference = 1) {
    .isRGOrStop(rgSet)
    normalize <- match.arg(normalize)

    if(normalize == "controls") {
        rgSet <- normalize.illumina.control(rgSet, reference = reference)
    }
    if(bg.correct) {
        rgSet <- bgcorrect.illumina(rgSet)
    }
    out <- preprocessRaw(rgSet)
    preprocess <- sprintf("Illumina, bg.correct = %s, normalize = %s, reference = %d",
                          bg.correct, normalize, reference)
    ## TODO:
    ## The manifest package version is currently not updated since 'packageVersion(getManifest(rgSet))' fails.
    ## packageVersion expects a string
    out@preprocessMethod <- c(rg.norm = preprocess,
                              minfi = as.character(packageVersion("minfi")),
                              manifest = as.character(packageVersion(.getManifestString(rgSet@annotation))))
    #packageVersion(getManifest(rgSet)))
    out
}

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
