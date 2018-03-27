# TODO: Refactor using similar/same strategy as in preprocessRaw()
.detectionP_matrix <- function(dimnames, Red, Green, controlIdx, TypeI.Red,
                               TypeI.Green, TypeII) {
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

.detectionP_DelayedMatrix <- function(dimnames, Red, Green, controlIdx,
                                      TypeI.Red, TypeI.Green, TypeII,
                                      BACKEND = getRealizationBackend()) {
    # Set up intermediate RealizationSink object of appropriate dimensions and
    # type
    # NOTE: This is ultimately coerced to the output DelayedMatrix objects,
    #       `detP`
    output_dim <- lengths(dimnames)
    detP_sink <- DelayedArray:::RealizationSink(dim = output_dim,
                                                dimnames = dimnames,
                                                type = "double")
    # We're going to walk over columns of `Red` and `Green` and write to columns
    # of `detP_sink`. This imposes some requirements:
    # 1a. Need to increase the block length so each block is made of at least
    #     one column.
    # 1b. Have to load columns of `Red` and `Green` into memory at the same time,
    #     so maximum block length is half what it would otherwise be if we only
    #     had to load one into memory.
    # 2. The grid over `detP_sink` must have the same dim as the grids over `Red`
    #    and `Green`.
    max_block_len <- 0.5 * max(nrow(Red) + nrow(Green),
                               DelayedArray:::get_max_block_length("double"))
    Red_grid <- defaultGrid(Red, max_block_len)
    Green_grid <- defaultGrid(Green, max_block_len)
    detP_sink_grid <- RegularArrayGrid(
        refdim = output_dim,
        spacings = c(output_dim[[1L]], output_dim[[2L]] / length(Red_grid)))
    # Sanity check ArrayGrid objects have the same dim
    stopifnot(dim(detP_sink_grid) == dim(Red_grid))
    # Loop over blocks of `Red` and `Green` and write to `detP_sink`
    # TODO: Adapted from DelayedArray::blockApply(); could blockApply() or
    #       blockMapply() be used directly?
    nblock <- length(Red_grid)
    bplapply(seq_len(nblock), function(b) {
        if (DelayedArray:::get_verbose_block_processing()) {
            message("Processing block ", b, "/", nblock, " ... ",
                    appendLF = FALSE)
        }
        # Set up viewports
        Red_viewport <- Red_grid[[b]]
        Green_viewport <- Green_grid[[b]]
        detP_sink_viewport <- detP_sink_grid[[b]]

        # Extract blocks of `Red` and `Green` as ordinary arrays
        Red_block <- DelayedArray:::extract_block(Red, Red_grid[[b]])
        if (!is.matrix(Red_block)) {
            Red_block <- as.matrix(Red_block)
        }
        attr(Red_block, "from_grid") <- Red_grid
        attr(Red_block, "block_id") <- b
        Green_block <- DelayedArray:::extract_block(Green, Green_viewport)
        if (!is.matrix(Green_block)) {
            Green_block <- as.matrix(Green_block)
        }
        attr(Green_block, "from_grid") <- Green_grid
        attr(Green_block, "block_id") <- b

        # Extract dimnames of block(s)
        cols <- seq(start(Red_viewport)[2L], end(Red_viewport)[2L])
        dimnames_block <- dimnames
        dimnames_block[[2L]] <- dimnames[[2L]][cols]

        # Construct blocks of `detP` as ordinary arrays and write to `detP_sink`
        detP_block <- .detectionP_matrix(
            dimnames = dimnames_block,
            Red = Red_block,
            Green = Green_block,
            controlIdx = controlIdx,
            TypeI.Red = TypeI.Red,
            TypeI.Green = TypeI.Green,
            TypeII = TypeII)
        write_block_to_sink(detP_block, detP_sink, detP_sink_viewport)
        detP_block <- NULL
        if (DelayedArray:::get_verbose_block_processing()) {
            message("OK")
        }
    }, BPPARAM = SerialParam())

    # Coerce intermediate RealizationSink object to output DelayedMatrix object
    as(detP_sink, "DelayedArray")
}

# TODO: In the release version of minfi, detectionP() always returns a matrix.
#       This version can return a matrix or a DelayedMatrix. Is this a good
#       idea and/or necessary?
detectionP <- function(rgSet, type = "m+u") {
    .isRGOrStop(rgSet)

    # Extract data to pass to low-level function that constructs `detP`
    locusNames <- getManifestInfo(rgSet, "locusNames")
    sampleNames <- sampleNames(rgSet)
    dimnames <- list(locusNames, sampleNames)
    controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")
    Red <- getRed(rgSet)
    Green <- getGreen(rgSet)
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    TypeII <- getProbeInfo(rgSet, type = "II")

    # Construct `detP`
    if (is(Red, "DelayedMatrix") || is(Green, "DelayedMatrix")) {
        detP <- .detectionP_DelayedMatrix(dimnames, Red, Green, controlIdx,
                                          TypeI.Red, TypeI.Green, TypeII)
    } else if (is.matrix(Red) && is.matrix(Green)) {
        detP <- .detectionP_matrix(dimnames, Red, Green, controlIdx, TypeI.Red,
                                   TypeI.Green, TypeII)
    } else {
        stop("'Red' and 'Green' assays must be 'matrix' or 'DelayedMatrix' objects")
    }

    detP
}
