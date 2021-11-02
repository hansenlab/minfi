# Internal generics ------------------------------------------------------------

# `...` are additional arguments passed to methods.
setGeneric(
    ".preprocessRawMeth",
    function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII, ...)
        standardGeneric(".preprocessRawMeth"),
    signature = c("Red", "Green"))

# `...` are additional arguments passed to methods.
setGeneric(
    ".preprocessRawUnmeth",
    function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII, ...)
        standardGeneric(".preprocessRawUnmeth"),
    signature = c("Red", "Green"))

# `...` are additional arguments passed to methods.
setGeneric(
    ".preprocessRaw",
    function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII, ...)
        standardGeneric(".preprocessRaw"),
    signature = c("Red", "Green"))

# Internal methods -------------------------------------------------------------

setMethod(
    ".preprocessRawMeth",
    c("matrix", "matrix"),
    function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII) {
        # Set up output matrix with appropriate dimension and type
        type <- .highestType(Red, Green)
        M <- matrix(.NA_type(type),
                    nrow = length(locusNames),
                    ncol = ncol(Red),
                    dimnames = list(locusNames, colnames(Red)))

        # Fill output matrix
        M[TypeI.Red$Name, ] <- Red[TypeI.Red$AddressB, ]
        M[TypeI.Green$Name, ] <- Green[TypeI.Green$AddressB, ]
        M[TypeII$Name, ] <- Green[TypeII$AddressA, ]

        # Return output matrix
        M
    }
)

setMethod(
    ".preprocessRawUnmeth",
    c("matrix", "matrix"),
    function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII, ...) {
        # Set up output matrix with appropriate dimension and type
        type <- .highestType(Red, Green)
        U <- matrix(.NA_type(type),
                    nrow = length(locusNames),
                    ncol = ncol(Red),
                    dimnames = list(locusNames, colnames(Red)))

        # Fill output matrix
        U[TypeI.Red$Name, ] <- Red[TypeI.Red$AddressA, ]
        U[TypeI.Green$Name, ] <- Green[TypeI.Green$AddressA, ]
        U[TypeII$Name, ] <- Red[TypeII$AddressA, ]

        # Return output matrix
        U
    }
)

setMethod(
    ".preprocessRaw",
    c("matrix", "matrix"),
    function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII, ...) {
        list(
            M = .preprocessRawMeth(
                Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII),
            U = .preprocessRawUnmeth(
                Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII))
    }
)

setMethod(
    ".preprocessRaw",
    c("DelayedMatrix", "DelayedMatrix"),
    function(Red, Green, locusNames, TypeI.Red, TypeI.Green, TypeII,
             BPREDO = list(), BPPARAM = SerialParam()) {
        # Set up intermediate RealizationSink objects of appropriate dimensions
        # and type
        # NOTE: These are ultimately coerced to the output DelayedMatrix
        #       objects, `M` and `U`
        # NOTE: Don't do `U_sink <- M_sink` or else these will reference the
        #       same object and clobber each other when written to!
        ans_type <- .highestType(Red, Green)
        M_sink <- DelayedArray::AutoRealizationSink(
            dim = c(length(locusNames), ncol(Red)),
            dimnames = list(locusNames, colnames(Red)),
            type = ans_type)
        on.exit(close(M_sink))
        U_sink <- DelayedArray::AutoRealizationSink(
            dim = c(length(locusNames), ncol(Red)),
            dimnames = list(locusNames, colnames(Red)),
            type = ans_type)
        on.exit(close(U_sink), add = TRUE)

        # Set up ArrayGrid instances over `Red` and `Green` as well as
        # "parallel" ArrayGrid instances over `M_sink` and `U_sink`.
        Red_grid <- colAutoGrid(Red)
        Green_grid <- colAutoGrid(Green)
        M_sink_grid <- RegularArrayGrid(
            refdim = dim(M_sink),
            spacings = c(nrow(M_sink), Red_grid@spacings[2]))
        U_sink_grid <- M_sink_grid
        # Sanity check ArrayGrid objects have the same dim
        stopifnot(dim(Red_grid) == dim(Green_grid),
                  dim(Red_grid) == dim(M_sink_grid))

        # Loop over blocks of `Red` and `Green` and write to `M_sink` and
        # `U_sink`
        blockMapplyWithRealization(
            FUN = .preprocessRaw,
            Red = Red,
            Green = Green,
            MoreArgs = list(
                locusNames = locusNames,
                TypeI.Red = TypeI.Red,
                TypeI.Green = TypeI.Green,
                TypeII = TypeII),
            sinks = list(M_sink, U_sink),
            dots_grids = list(Red_grid, Green_grid),
            sinks_grids = list(M_sink_grid, U_sink_grid),
            BPREDO = BPREDO,
            BPPARAM = BPPARAM)

        # Return as DelayedMatrix objects
        M <- as(M_sink, "DelayedArray")
        U <- as(U_sink, "DelayedArray")

        list(M = M, U = U)
    }
)

# Exported functions -----------------------------------------------------------

# TODO: Document: because we simultaneously walk over column-blocks of
#       `Red` and `Green`, the number of elements loaded into memory is
#       doubled. If running into memory issues, try halving
#       getOption("DelayedArray.block.size")
preprocessRaw <- function(rgSet) {
    .isRGOrStop(rgSet)

    # Extract data to pass to low-level functions that construct `M` and `U`
    locusNames <- getManifestInfo(rgSet, "locusNames")
    Red <- getRed(rgSet)
    Green <- getGreen(rgSet)
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    TypeII <- getProbeInfo(rgSet, type = "II")

    # Construct `M` and `U`
    M_and_U <- .preprocessRaw(
        Red = Red,
        Green = Green,
        locusNames = locusNames,
        TypeI.Red = TypeI.Red,
        TypeI.Green = TypeI.Green,
        TypeII = TypeII)
    M <- M_and_U[["M"]]
    U <- M_and_U[["U"]]

    # Construct MethylSet
    out <- MethylSet(
        Meth = M,
        Unmeth = U,
        colData = colData(rgSet),
        annotation = annotation(rgSet),
        metadata = metadata(rgSet))

    # TODO: The manifest package version is currently not updated since
    #       `packageVersion(getManifest(rgSet))` fails. `packageVersion()
    #       expects a string
    out@preprocessMethod <- c(
        rg.norm = "Raw (no normalization or bg correction)",
        minfi = as.character(packageVersion("minfi")),
        manifest = as.character(
            packageVersion(.getManifestString(rgSet@annotation))))
    #packageVersion(getManifest(rgSet)))
    out
}
