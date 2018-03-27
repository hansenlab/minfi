.preprocessRawMeth <- function(dimnames, Red, Green, TypeI.Red, TypeI.Green,
                               TypeII) {
    # Set up output matrix with appropriate dimension and type
    dim <- lengths(dimnames)
    type <- .highestType(Red, Green)
    M <- matrix(.NA_type(type),
                nrow = dim[[1L]],
                ncol = dim[[2L]],
                dimnames = dimnames)

    # Fill output matrix
    M[TypeI.Red$Name, ] <- Red[TypeI.Red$AddressB, ]
    M[TypeI.Green$Name, ] <- Green[TypeI.Green$AddressB, ]
    M[TypeII$Name, ] <- Green[TypeII$AddressA, ]

    # Return output matrix
    M
}

.preprocessRawUnmeth <- function(dimnames, Red, Green, TypeI.Red, TypeI.Green,
                                 TypeII) {
    # Set up output matrix with appropriate dimension and type
    dim <- lengths(dimnames)
    type <- .highestType(Red, Green)
    U <- matrix(.NA_type(type),
                nrow = dim[[1L]],
                ncol = dim[[2L]],
                dimnames = dimnames)

    # Fill output matrix
    U[TypeI.Red$Name, ] <- Red[TypeI.Red$AddressA, ]
    U[TypeI.Green$Name, ] <- Green[TypeI.Green$AddressA, ]
    U[TypeII$Name, ] <- Red[TypeII$AddressA, ]

    # Return output matrix
    U
}

.preprocessRaw_DelayedMatrix <- function(dimnames, Red, Green, TypeI.Red,
                                         TypeI.Green, TypeII) {
    # Set up intermediate RealizationSink objects of appropriate dimensions and
    # type
    # NOTE: These are ultimately coerced to the output DelayedMatrix objects,
    #       `M` and `U`
    output_dim <- lengths(dimnames)
    type <- .highestType(Red, Green)
    M_sink <- DelayedArray:::RealizationSink(dim = output_dim,
                                             dimnames = dimnames,
                                             type = type)
    # NOTE: Don't do `U_sink <- M_sink` or else these will reference the same
    #       object and clobber each other when written to!
    U_sink <- DelayedArray:::RealizationSink(dim = output_dim,
                                             dimnames = dimnames,
                                             type = type)

    # Set up ArrayGrid instances over `Red`, `Green`, `M_sink`, and `U_sink`.
    # We're going to walk over columns of `Red` and `Green` and write to columns
    # of `M_sink` and `U_sink`. This imposes some requirements:
    # 1a. Need to increase the block length so each block is made of at least
    #     one column.
    # 1b. Have to load columns of `Red` and `Green` into memory at the same time,
    #     so maximum block length is half what it would otherwise be if we only
    #     had to load one into memory.
    # 2. The grids over `M_sink` and `U_sink` must have the same dim as the grids
    #    over `Red` and `Green`.
    max_block_len <- 0.5 * max(nrow(Red) + nrow(Green),
                               DelayedArray:::get_max_block_length(type))
    Red_grid <- defaultGrid(Red, max_block_len)
    Green_grid <- defaultGrid(Green, max_block_len)
    M_sink_grid <- RegularArrayGrid(
        refdim = output_dim,
        spacings = c(output_dim[[1L]], output_dim[[2L]] / length(Red_grid)))
    U_sink_grid <- M_sink_grid
    # Sanity check ArrayGrid objects have the same dim
    stopifnot(dim(M_sink_grid) == dim(Red_grid))

    # Loop over blocks of `Red` and `Green` and write to `M_sink` and `U_sink`
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
        M_sink_viewport <- M_sink_grid[[b]]
        U_sink_viewport <- U_sink_grid[[b]]

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

        # Construct blocks of `M`  and `U` as ordinary arrays and write
        # to `M_sink` and `U_sink`, respectively.
        M_block <- .preprocessRawMeth(
            dimnames = dimnames_block,
            Red = Red_block,
            Green = Green_block,
            TypeI.Red = TypeI.Red,
            TypeI.Green = TypeI.Green,
            TypeII = TypeII)
        write_block_to_sink(M_block, M_sink, M_sink_viewport)
        M_block <- NULL
        U_block <- .preprocessRawUnmeth(
            dimnames = dimnames_block,
            Red = Red_block,
            Green = Green_block,
            TypeI.Red = TypeI.Red,
            TypeI.Green = TypeI.Green,
            TypeII = TypeII)
        write_block_to_sink(U_block, U_sink, U_sink_viewport)
        U_block <- NULL
        if (DelayedArray:::get_verbose_block_processing()) {
            message("OK")
        }
    }, BPPARAM = SerialParam())

    # Coerce intermediate RealizationSink objects to output DelayedMatrix objects
    M <- as(M_sink, "DelayedArray")
    U <- as(U_sink, "DelayedArray")

    list(M = M, U = U)
}

# TODO: Add `BACKEND` argument so that output can have different backend to input?
#       Or should that be driven by setRealizationBackend()?
# TODO: Add BPPARAM (and other BiocParallel args, e.g., BPREDO?)?
# TODO: Instead of (or in addition to) relying on option("DelayedArray.block.size"),
#       perhaps offer argument to specify how many samples should be processed per
#       "batch".
preprocessRaw <- function(rgSet) {
    .isRGOrStop(rgSet)

    # Extract data to pass to low-level functions that construct `M` and `U`
    locusNames <- getManifestInfo(rgSet, "locusNames")
    sampleNames <- sampleNames(rgSet)
    dimnames <- list(locusNames, sampleNames)
    Red <- getRed(rgSet)
    Green <- getGreen(rgSet)
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    TypeII <- getProbeInfo(rgSet, type = "II")

    # Construct `M` and `U`
    if (is(Red, "DelayedMatrix") || is(Green, "DelayedMatrix")) {
        M_and_U <- .preprocessRaw_DelayedMatrix(dimnames, Red, Green, TypeI.Red,
                                                TypeI.Green, TypeII)
        M <- M_and_U[["M"]]
        U <- M_and_U[["U"]]
    } else if (is.matrix(Red) && is.matrix(Green)) {
        M <- .preprocessRawMeth(dimnames, Red, Green, TypeI.Red, TypeI.Green, TypeII)
        U <- .preprocessRawUnmeth(dimnames, Red, Green, TypeI.Red, TypeI.Green, TypeII)
    } else {
        stop("'Red' and 'Green' assays must be 'matrix' or 'DelayedMatrix' objects")
    }

    # Construct MethylSet
    out <- MethylSet(Meth = M, Unmeth = U, colData = colData(rgSet),
                     annotation = annotation(rgSet), metadata = metadata(rgSet))
    ## TODO:
    ## The manifest package version is currently not updated since `packageVersion(getManifest(rgSet))` fails.
    ## packageVersion expects a string
    out@preprocessMethod <- c(
        rg.norm = "Raw (no normalization or bg correction)",
        minfi = as.character(packageVersion("minfi")),
        manifest = as.character(packageVersion(.getManifestString(rgSet@annotation))))
    #packageVersion(getManifest(rgSet)))
    out
}
