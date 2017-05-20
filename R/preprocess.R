# TODO: Add mc.cores to parallelise across samples?
preprocessRaw <- function(rgSet) {
    .isRGOrStop(rgSet)
    if (.isHDF5ArrayBacked(rgSet)) {
        BACKEND <- "HDF5Array"
    } else {
        BACKEND <- NULL
    }
    locusNames <- getManifestInfo(rgSet, "locusNames")
    # HDF5: Idea is to process sample-by-sample (column-by-column). Each
    #       iteration loads the red and green data into memory as a matrix,
    #       then subsets  by rows, then writes back to disk. This is an
    #       efficient strategy because each sample's data is contiguous in the
    #       HDF5 file. If instead we first did rowwise subsetting, then
    #       reading the data from the HDF5 file incurs a large overhead due to
    #       non-contiguous reads.
    TypeII <- getProbeInfo(rgSet, type = "II")
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    M_and_U <- lapply(sampleNames(rgSet), function(this_sample) {
        M <- matrix(NA_real_, ncol = 1, nrow = length(locusNames),
                    dimnames = list(locusNames, this_sample))
        U <- M
        green <- as.matrix(getGreen(rgSet)[, this_sample])
        red <- as.matrix(getRed(rgSet)[, this_sample])
        M[TypeII.Name, ] <- green[TypeII$AddressA, ]
        U[TypeII.Name, ] <- red[TypeII$AddressA, ]
        M[TypeI.Red$Name,] <- red[TypeI.Red$AddressB, ]
        M[TypeI.Green$Name,] <- green[TypeI.Green$AddressB, ]
        U[TypeI.Red$Name,] <- red[TypeI.Red$AddressA, ]
        U[TypeI.Green$Name,] <- green[TypeI.Green$AddressA, ]
        # NOTE: Don't store dimnames in DelayedArray (these get added back
        #       below once samples are cbind()-ed)
        M <- realize(unname(M), BACKEND = BACKEND)
        U <- realize(unname(U), BACKEND = BACKEND)
        list(M = M, U = U)
    })
    M <- do.call(cbind, lapply(M_and_U, "[[", "M"))
    U <- do.call(cbind, lapply(M_and_U, "[[", "U"))
    # NOTE: Add back the dimnames
    dimnames(M) <- list(locusNames, sampleNames(rgSet))
    dimnames(U) <- list(locusNames, sampleNames(rgSet))
    out <- MethylSet(Meth = M,
                     Unmeth = U,
                     colData = colData(rgSet),
                     annotation = annotation(rgSet),
                     metadata = metadata(rgSet))
    ## TODO:
    ## The manifest package version is currently not updated since 'packageVersion(getManifest(rgSet))' fails.
    ## packageVersion expects a string
    out@preprocessMethod <- c(
        rg.norm = "Raw (no normalization or bg correction)",
        minfi = as.character(packageVersion("minfi")),
        manifest = as.character(packageVersion(
            .getManifestString(rgSet@annotation))))
    #packageVersion(getManifest(rgSet)))
    out
}

# HDF5: Not yet supported
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
    # TODO: Need a sweep() for DelayedArray
    Green <- sweep(Green, 2, FUN = "*", Green.factor)
    Red <- sweep(Red, 2, FUN = "*", Red.factor)
    assay(rgSet, "Green") <- Green
    assay(rgSet, "Red") <- Red
    rgSet
}

# HDF5: Not yet supported
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
    # NOTE: sort() doesn't work on DelayedArray, so need to realize each
    #       column as a matrix
    Green.bg <- apply(Green[NegControls, , drop = FALSE], 2, function(xx) {
        sort(as.matrix(xx))[31]
        })
    Red.bg <- apply(Red[NegControls, , drop = FALSE], 2, function(xx) {
        sort(as.matrix(xx))[31]
    })
    # TODO: Need a sweep() for DelayedArray
    Green <- pmax(sweep(Green, 2, Green.bg), 0)
    Red <- pmax(sweep(Red, 2, Red.bg), 0)
    assay(rgSet, "Green") <- Green
    assay(rgSet, "Red") <- Red
    rgSet
}


# HDF5: Not yet supported
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

# TODO: Add mc.cores to parallelise across samples?
detectionP <- function(rgSet, type = "m+u") {
    .isRGOrStop(rgSet)
    if (.isHDF5ArrayBacked(rgSet)) {
        BACKEND <- "HDF5Array"
    } else {
        BACKEND <- NULL
    }
    locusNames <- getManifestInfo(rgSet, "locusNames")
    controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")
    TypeII <- getProbeInfo(rgSet, type = "II")
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    # HDF5: Idea is to process sample-by-sample (column-by-column). Each
    #       iteration loads the red and green data into memory as a matrix,
    #       then subsets  by rows, then writes back to disk. This is an
    #       efficient strategy because each sample's data is contiguous in the
    #       HDF5 file. If instead we first did rowwise subsetting, then
    #       reading the data from the HDF5 file incurs a large overhead due to
    #       non-contiguous reads.
    detP <- do.call(cbind, lapply(sampleNames(rgSet), function(this_sample) {
        detP <- matrix(NA_real_, ncol = 1, nrow = length(locusNames),
                       dimnames = list(locusNames, this_sample))

        r <- as.matrix(getRed(rgSet)[, j, drop = FALSE])
        rBg <- r[controlIdx, , drop = FALSE]
        # NOTE: matrixStats::colMedians() is overkill since this is a
        #       one-column matrix, but it's still faster than median()
        rMu <- matrixStats::colMedians(rBg)
        rSd <- matrixStats::colMads(rBg)

        g <- as.matrix(getGreen(rgSet)[, j, drop = FALSE])
        gBg <- g[controlIdx, , drop = FALSE]
        gMu <- matrixStats::colMedians(gBg)
        gSd <- matrixStats::colMads(gBg)

        ## Type I Red
        intensity <- r[TypeI.Red$AddressA, j] + r[TypeI.Red$AddressB, j]
        detP[TypeI.Red$Name, ] <- 1 - pnorm(intensity,
                                            mean = rMu * 2,
                                            sd = rSd * 2)
        ## Type I Green
        intensity <- g[TypeI.Green$AddressA, ] + g[TypeI.Green$AddressB, ]
        detP[TypeI.Green$Name, ] <- 1 - pnorm(intensity,
                                              mean = gMu * 2,
                                              sd = gSd * 2)
        ## Type II
        intensity <- r[TypeII$AddressA, ] + g[TypeII$AddressA, ]
        detP[TypeII$Name, ] <- 1 - pnorm(intensity,
                                         mean = rMu + gMu,
                                         sd = rSd + gSd)

        # NOTE: Don't store dimnames in DelayedArray (these get added back
        #       below once samples are cbind()-ed)
        realize(unname(detP), BACKEND = BACKEND)
    }))
    dimnames(detP) <- list(locusNames, sampleNames(rgSet))
    detP
}
