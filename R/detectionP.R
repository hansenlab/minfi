# Internal generics ------------------------------------------------------------

# NOTE: `...` are additional arguments passed to methods.
setGeneric(
    ".detectionP",
    function(Red, Green, locusNames, controlIdx, TypeI.Red, TypeI.Green, TypeII,
             ...)
        standardGeneric(".detectionP"),
    signature = c("Red", "Green"))

# Internal methods -------------------------------------------------------------

# NOTE: `...` is ignored
setMethod(
    ".detectionP",
    c("matrix", "matrix"),
    function(Red, Green, locusNames, controlIdx, TypeI.Red, TypeI.Green, TypeII,
             ...) {
        # Set up output matrix with appropriate dimensions and type
        detP <- matrix(NA_real_,
                       nrow = length(locusNames),
                       ncol = ncol(Red),
                       dimnames = list(locusNames, colnames(Red)))

        # Compute summary statistics needed for calculations
        rBg <- Red[controlIdx, , drop = FALSE]
        rMu <- colMedians(rBg)
        rSd <- colMads(rBg)
        gBg <- Green[controlIdx, , drop = FALSE]
        gMu <- colMedians(gBg)
        gSd <- colMads(gBg)

        # Fill output matrix
        for (j in seq_len(ncol(detP))) {
            # Type I Red
            intensity <- Red[TypeI.Red$AddressA, j] + Red[TypeI.Red$AddressB, j]
            detP[TypeI.Red$Name, j] <- pnorm(
                q = intensity,
                mean = 2 * rMu[j],
                sd = 2 * rSd[j],
                lower.tail = FALSE)
            # Type I Green
            intensity <- Green[TypeI.Green$AddressA, j] +
                Green[TypeI.Green$AddressB, j]
            detP[TypeI.Green$Name, j] <- pnorm(
                q = intensity,
                mean = 2 * gMu[j],
                sd = 2 * gSd[j],
                lower.tail = FALSE)
            # Type II
            intensity <- Red[TypeII$AddressA, j] + Green[TypeII$AddressA, j]
            detP[TypeII$Name, j] <- pnorm(
                q = intensity,
                mean = rMu[j] + gMu[j],
                sd = rSd[j] + gSd[j],
                lower.tail = FALSE)
        }

        # Return output matrix
        detP
    }
)

setMethod(
    ".detectionP",
    c("DelayedMatrix", "DelayedMatrix"),
    function(Red, Green, locusNames, controlIdx, TypeI.Red, TypeI.Green, TypeII,
             BPREDO = list(), BPPARAM = SerialParam()) {
        # Set up intermediate RealizationSink object of appropriate dimensions
        # and type
        # NOTE: This is ultimately coerced to the output DelayedMatrix objects,
        #       `detP`
        detP_sink <- DelayedArray::AutoRealizationSink(
            dim = c(length(locusNames), ncol(Red)),
            dimnames = list(locusNames, colnames(Red)),
            type = "double")
        on.exit(close(detP_sink))

        # Set up ArrayGrid instances over `Red` and `Green` as well as
        # "parallel" ArrayGrid instance over `detP_sink`.
        Red_grid <- colAutoGrid(Red)
        Green_grid <- colAutoGrid(Green)
        detP_sink_grid <- RegularArrayGrid(
            refdim = dim(detP_sink),
            spacings = c(nrow(detP_sink), ncol(detP_sink) / length(Red_grid)))
        # Sanity check ArrayGrid objects have the same dim
        stopifnot(dim(Red_grid) == dim(Green_grid),
                  dim(Red_grid) == dim(detP_sink_grid))

        # Loop over blocks of `Red` and `Green` and write to `detP_sink`
        blockMapplyWithRealization(
            FUN = .detectionP,
            Red = Red,
            Green = Green,
            MoreArgs = list(
                locusNames = locusNames,
                controlIdx = controlIdx,
                TypeI.Red = TypeI.Red,
                TypeI.Green = TypeI.Green,
                TypeII = TypeII),
            sinks = list(detP_sink),
            dots_grids = list(Red_grid, Green_grid),
            sinks_grids = list(detP_sink_grid),
            BPREDO = BPREDO,
            BPPARAM = BPPARAM)

        # Return as DelayedMatrix
        as(detP_sink, "DelayedArray")
    }
)

# Exported functions -----------------------------------------------------------

# TODO: Document: because we simultaneously walk over column-blocks of
#       `Red` and `Green`, the number of elements loaded into memory is
#       doubled. If running into memory issues, try halving
#       getOption("DelayedArray.block.size")
detectionP <- function(rgSet, type = "m+u") {
    .isRGOrStop(rgSet)

    # Extract data to pass to low-level function that constructs `detP`
    locusNames <- getManifestInfo(rgSet, "locusNames")
    controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")
    Red <- getRed(rgSet)
    Green <- getGreen(rgSet)
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    TypeII <- getProbeInfo(rgSet, type = "II")

    # Construct `detP`
    detP <- .detectionP(
        Red = Red,
        Green = Green,
        locusNames = locusNames,
        controlIdx = controlIdx,
        TypeI.Red = TypeI.Red,
        TypeI.Green = TypeI.Green,
        TypeII = TypeII)

    detP
}
