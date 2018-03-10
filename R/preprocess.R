.preProcessRaw_matrix <- function(red, green, locusNames, sampleNames, TypeI.Red,
                                  TypeI.Green, TypeII) {
    M <- matrix(NA_real_,
                ncol = ncol(rgSet),
                nrow = length(locusNames),
                dimnames = list(locusNames, sampleNames))
    U <- M

    M[TypeII$Name, ] <- green[TypeII$AddressA, ]
    U[TypeII$Name, ] <- red[TypeII$AddressA, ]

    M[TypeI.Red$Name, ] <- red[TypeI.Red$AddressB, ]
    M[TypeI.Green$Name, ] <- green[TypeI.Green$AddressB, ]
    U[TypeI.Red$Name, ] <- red[TypeI.Red$AddressA, ]
    U[TypeI.Green$Name, ] <- green[TypeI.Green$AddressA, ]

    list(M = M, U = U)
}

.preProcessRaw_DelayedArray <- function(rgSet, locusNames, TypeI.Red, TypeI.Green,
                                        TypeII, BACKEND) {
    nrow <- length(locusNames)
    sampleNames <- sampleNames(rgSet)
    ncol <- length(sampleNames)
    red <- getRed(rgSet)
    green <- getGreen(rgSet)
    # TODO: Take the 'higher' type of red and green
    type <- DelayedArray:::type(red)

    # NOTE:  We're going to walk along the columns so need to increase the block
    #        length so each block is made of at least one column.
    max_block_len <- max(
        DelayedArray:::get_max_block_length(type), nrow)
    red_grid <- defaultGrid(red, max_block_len)
    green_grid <- defaultGrid(green, max_block_len)

    ### ------------------------------------------------------------------------------
    ### M
    ###

    # TODO: It feels like 'M' shouldn't be necessary
    rle <- Rle(vector(type, 1L), prod(nrow, ncol))
    M <- RleArray(rle,
                  dim = c(nrow, ncol),
                  dimnames = list(locusNames, sampleNames))
    # TODO: What is the proper way to set up the realization sink to match the
    #       BACKEND? Within this function, may need to override the global value
    M_sink <- DelayedArray:::RealizationSink(
        dim = c(nrow, ncol),
        dimnames = list(locusNames, sampleNames),
        type = type)
    M_grid <- RegularArrayGrid(dim(M), c(nrow(M), ncol(M) / ncol(red_grid)))
    stopifnot(length(M_grid) == length(red_grid))
    nblock <- length(M_grid)

    M_blocks <- bplapply(
        seq_len(nblock),
        function(b, red_grid, green_grid, M_grid, red, green, M, TypeI.Red,
                 TypeI.Green, TypeII) {
            if (DelayedArray:::get_verbose_block_processing()) {
                message("Processing block ", b, "/", nblock, " ... ",
                        appendLF = FALSE)
            }
            red_viewport <- red_grid[[b]]
            green_viewport <- green_grid[[b]]
            M_viewport <- M_grid[[b]]
            red_block <- DelayedArray:::extract_block(red, red_viewport)
            green_block <- DelayedArray:::extract_block(green, green_viewport)
            M_block <- DelayedArray:::extract_block(M, M_viewport)
            if (!is.array(red_block)) {
                red_block <- DelayedArray:::.as_array_or_matrix(red_block)
            }
            attr(red_block, "from_grid") <- red_grid
            attr(red_block, "block_id") <- b
            if (!is.array(green_block)) {
                green_block <- DelayedArray:::.as_array_or_matrix(green_block)
            }
            attr(green_block, "from_grid") <- green_grid
            attr(green_block, "block_id") <- b
            if (!is.array(M_block)) {
                M_block <- DelayedArray:::.as_array_or_matrix(M_block)
            }
            attr(M_block, "from_grid") <- M_grid
            attr(M_block, "block_id") <- b

            # Update M
            M_block[TypeII$Name, ] <- green_block[TypeII$AddressA, ]
            M_block[TypeI.Red$Name, ] <- red_block[TypeI.Red$AddressB, ]
            M_block[TypeI.Green$Name, ] <- green_block[TypeI.Green$AddressB, ]

            if (!is.null(M_sink)) {
                write_block_to_sink(M_block, M_sink, M_viewport)
                M_block <- NULL
            }
            if (DelayedArray:::get_verbose_block_processing()) {
                message("OK")
            }
            M_block
        }, red_grid, green_grid, M_grid, red, green, M, TypeI.Red, TypeI.Green,
        TypeII)

    ### ------------------------------------------------------------------------------
    ### U
    ###U
    # TODO: It feels like 'M' shouldn't be necessary
    rle <- Rle(vector(type, 1L), prod(nrow, ncol))
    U <- RleArray(rle,
                  dim = c(nrow, ncol),
                  dimnames = list(locusNames, sampleNames))
    # TODO: What is the proper way to set up the realization sink to match the
    #       BACKEND? Within this function, may need to override the global value
    U_sink <- DelayedArray:::RealizationSink(
        dim = c(nrow, ncol),
        dimnames = list(locusNames, sampleNames),
        type = type)
    U_grid <- RegularArrayGrid(dim(U), c(nrow(U), ncol(U) / ncol(red_grid)))
    stopifnot(length(U_grid) == length(red_grid))
    nblock <- length(U_grid)
    # TODO: This successfully writes to the sink if the loop is run
    #       'manually' and when bpparam is a SerialParam, but not when it is
    #       a MultiCoreParam; I think it's because some variables aren't
    #       being exported to the workers
    U_blocks <- bplapply(
        seq_len(nblock),
        function(b, red_grid, green_grid, U_grid, red, green, M, TypeI.Red,
                 TypeI.Green,
                 TypeII) {
            if (DelayedArray:::get_verbose_block_processing()) {
                message("Processing block ", b, "/", nblock, " ... ",
                        appendLF = FALSE)
            }
            red_viewport <- red_grid[[b]]
            green_viewport <- green_grid[[b]]
            U_viewport <- U_grid[[b]]
            red_block <- DelayedArray:::extract_block(red, red_viewport)
            green_block <- DelayedArray:::extract_block(green, green_viewport)
            U_block <- DelayedArray:::extract_block(U, U_viewport)
            if (!is.array(red_block)) {
                red_block <- DelayedArray:::.as_array_or_matrix(red_block)
            }
            attr(red_block, "from_grid") <- red_grid
            attr(red_block, "block_id") <- b
            if (!is.array(green_block)) {
                green_block <- DelayedArray:::.as_array_or_matrix(green_block)
            }
            attr(green_block, "from_grid") <- green_grid
            attr(green_block, "block_id") <- b
            if (!is.array(U_block)) {
                U_block <- DelayedArray:::.as_array_or_matrix(U_block)
            }
            attr(U_block, "from_grid") <- U_grid
            attr(U_block, "block_id") <- b

            # Update U
            U_block[TypeII$Name, ] <- red_block[TypeII$AddressA, ]
            U_block[TypeI.Red$Name, ] <- red_block[TypeI.Red$AddressA, ]
            U_block[TypeI.Green$Name, ] <- green_block[TypeI.Green$AddressA, ]

            if (!is.null(U_sink)) {
                write_block_to_sink(U_block, U_sink, U_viewport)
                U_block <- NULL
            }
            if (DelayedArray:::get_verbose_block_processing()) {
                message("OK")
            }
            U_block
        }, red_grid, green_grid, U_grid, red, green, U, TypeI.Red, TypeI.Green,
        TypeII)

    if (is.null(BACKEND)) {
        M <- do.call(cbind, M_blocks)
        U <- do.call(cbind, U_blocks)
    } else {
        M <- as(M_sink, "DelayedArray")
        U <- as(U_sink, "DelayedArray")
    }
    list(M = M, U = U)
}

preprocessRaw <- function(rgSet, BACKEND = getRealizationBackend()) {
    .isRGOrStop(rgSet)

    # Extract data to pass to low-level preprocessing function
    red <- getRed(rgSet)
    green <- getGreen(rgSet)
    locusNames <- getManifestInfo(rgSet, "locusNames")
    sampleNames <- sampleNames(rgSet)
    TypeII <- getProbeInfo(rgSet, type = "II")
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")

    # Fill M and U matrices using probe names and info
    if (.isBackedByDelayedArray(rgSet)) {
        val <- .preProcessRaw_DelayedArray(rgSet, locusNames, TypeI.Red,
                                           TypeI.Green, TypeII, BACKEND)
    } else {
        val <- .preProcessRaw_matrix(red, green, locusNames, sampleNames, TypeI.Red,
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
    Green.bg <- apply(Green[NegControls, , drop = FALSE], 2, function(xx) sort(xx)[31])
    Red.bg <- apply(Red[NegControls, , drop = FALSE], 2, function(xx) sort(xx)[31])
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


detectionP <- function(rgSet, type = "m+u") {
    .isRGOrStop(rgSet)
    locusNames <- getManifestInfo(rgSet, "locusNames")
    detP <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
                   dimnames = list(locusNames, colnames(rgSet)))

    controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")
    r <- getRed(rgSet)
    rBg <- r[controlIdx,,drop=FALSE]
    rMu <- matrixStats::colMedians(rBg)
    rSd <- matrixStats::colMads(rBg)

    g <- getGreen(rgSet)
    gBg <- g[controlIdx,,drop=FALSE]
    gMu <- matrixStats::colMedians(gBg)
    gSd <- matrixStats::colMads(gBg)

    TypeII <- getProbeInfo(rgSet, type = "II")
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    for (i in 1:ncol(rgSet)) {
        ## Type I Red
        intensity <- r[TypeI.Red$AddressA, i] + r[TypeI.Red$AddressB, i]
        detP[TypeI.Red$Name, i] <- 1-pnorm(intensity, mean=rMu[i]*2, sd=rSd[i]*2)
        ## Type I Green
        intensity <- g[TypeI.Green$AddressA, i] + g[TypeI.Green$AddressB, i]
        detP[TypeI.Green$Name, i] <- 1-pnorm(intensity, mean=gMu[i]*2, sd=gSd[i]*2)
        ## Type II
        intensity <- r[TypeII$AddressA, i] + g[TypeII$AddressA, i]
        detP[TypeII$Name, i] <- 1-pnorm(intensity, mean=rMu[i]+gMu[i], sd=rSd[i]+gSd[i])
    }
    detP
}
