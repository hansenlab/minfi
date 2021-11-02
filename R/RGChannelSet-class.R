# Exported classes -------------------------------------------------------------

setClass(
    "RGChannelSet",
    representation(annotation = "character"),
    contains = "SummarizedExperiment"
)

# Validity methods -------------------------------------------------------------

setValidity("RGChannelSet", function(object) {
    msg <- validMsg(NULL, .checkAssayNames(object, c("Red", "Green")))
    if (is.null(msg)) TRUE else msg
})

# Internal generics ------------------------------------------------------------

# `x` and `probe_info` are either `Red` and `IGrn` or `Green` and `IRed`,
# respectively.
# `...` are additional arguments passed to methods.
setGeneric(
    ".getOOB",
    function(x, probe_info, ...) standardGeneric(".getOOB"),
    signature = c("x")
)

# Internal methods -------------------------------------------------------------

setMethod(".getOOB", "matrix", function(x, probe_info, ...) {
    rbind(
        x[probe_info$AddressA, , drop = FALSE],
        x[probe_info$AddressB, , drop = FALSE])
})

setMethod(
    ".getOOB",
    "DelayedMatrix",
    function(x, probe_info, BPREDO = list(), BPPARAM = SerialParam()) {
        # Set up intermediate RealizationSink object of appropriate dimensions
        # and type.
        # NOTE: This is ultimately coerced to the output DelayedMatrix object.
        ans_type <- type(x)
        sink <- DelayedArray::AutoRealizationSink(
            dim = c(2L * nrow(probe_info), ncol(x)),
            dimnames = list(
                c(probe_info$AddressA, probe_info$AddressB), colnames(x)),
            type = ans_type)
        on.exit(close(sink))

        # Set up ArrayGrid instances over `x` as well as "parallel" ArrayGrid
        # instance over `sink`.
        x_grid <- colAutoGrid(x)
        sink_grid <- RegularArrayGrid(
            refdim = dim(sink),
            spacings = c(nrow(sink), ncol(sink) / length(x_grid)))
        # Sanity check ArrayGrid objects have the same dim
        stopifnot(dim(x_grid) == dim(sink_grid))

        # Loop over blocks of `Green` and write to `Grn_sink`
        blockApplyWithRealization(
            FUN = .getOOB,
            x = x,
            probe_info = probe_info,
            sink = sink,
            x_grid = x_grid,
            sink_grid = sink_grid,
            BPREDO = BPREDO,
            BPPARAM = BPPARAM)

        # Return as DelayedMatrix
        as(sink, "DelayedArray")
    }
)

# Exported functions -----------------------------------------------------------

RGChannelSet <- function(Green = new("matrix"), Red = new("matrix"),
                         annotation = "", ...) {
    ## Check rownames, colnames
    assays <- SimpleList(Green = Green, Red = Red)
    new("RGChannelSet",
        SummarizedExperiment(assays = assays, ...),
        annotation = annotation
    )
}

getRed <- function(object) {
    .isRGOrStop(object)
    assay(object, "Red")
}

getGreen <- function(object) {
    .isRGOrStop(object)
    assay(object, "Green")
}

getOOB <- function(object) {
    # Check input
    .isRGOrStop(object)

    # Extract data to pass to low-level functions that construct `Red` and `Grn`
    Red <- getRed(object)
    Green <- getGreen(object)
    IRed <- getProbeInfo(object, type = "I-Red")
    IGreen <- getProbeInfo(object, type = "I-Green")

    # Extract OOB probes
    list(Grn = .getOOB(Green, IRed), Red = .getOOB(Red, IGreen))
}

getSnpBeta <- function(object) {
    # Check input
    .isRGOrStop(object)

    # Extract and processSNP probes
    snpProbesI <- getProbeInfo(object, type = "SnpI")
    snpProbesII <- getProbeInfo(object, type = "SnpII")
    snpProbesI.Green <- snpProbesI[snpProbesI$Color == "Grn", , drop = FALSE]
    snpProbesI.Red <- snpProbesI[snpProbesI$Color == "Red", , drop = FALSE]

    # Extract Red and Green channels
    Green <- getGreen(object)
    Red <- getRed(object)

    # Construct M and U
    M.II <- Green[snpProbesII$AddressA, , drop = FALSE]
    U.II <- Red[snpProbesII$AddressA, , drop = FALSE]
    M.I.Red <- Red[snpProbesI.Red$AddressB, , drop = FALSE]
    U.I.Red <- Red[snpProbesI.Red$AddressA, , drop = FALSE]
    M.I.Green <- Green[snpProbesI.Green$AddressB, , drop = FALSE]
    U.I.Green <- Green[snpProbesI.Green$AddressA, , drop = FALSE]
    M <- rbind(M.II, M.I.Red, M.I.Green)
    U <- rbind(U.II, U.I.Red, U.I.Green)
    rownames(M) <- rownames(U) <- c(
        snpProbesII$Name,
        snpProbesI.Red$Name,
        snpProbesI.Green$Name)

    # Compute beta
    M / (U + M + 100)
}

subsetByLoci <- function(rgSet, includeLoci = NULL, excludeLoci = NULL,
                         keepControls = TRUE, keepSnps = TRUE) {
    # Check input
    .isRGOrStop(rgSet)

    # Get probe info
    TypeII <- getProbeInfo(rgSet, type = "II")
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")

    # Select probes to retain
    if (!is.null(includeLoci)) {
        TypeII <- TypeII[TypeII$Name %in% includeLoci, ]
        TypeI.Red <- TypeI.Red[TypeI.Red$Name %in% includeLoci, ]
        TypeI.Green <- TypeI.Green[TypeI.Green$Name %in% includeLoci, ]
    }
    if (!is.null(excludeLoci)) {
        TypeII <- TypeII[!TypeII$Name %in% excludeLoci, ]
        TypeI.Red <- TypeI.Red[!TypeI.Red$Name %in% excludeLoci, ]
        TypeI.Green <- TypeI.Green[!TypeI.Green$Name %in% excludeLoci, ]
    }
    addresses <- c(
        TypeII$AddressA,
        TypeI.Red$AddressA, TypeI.Red$AddressB,
        TypeI.Green$AddressA, TypeI.Green$AddressB)
    if (keepControls) {
        addresses <- c(addresses, getProbeInfo(rgSet, type = "Control")$Address)
    }
    if (keepSnps) {
        addresses <- c(
            addresses,
            getProbeInfo(rgSet, type = "SnpI")$AddressA,
            getProbeInfo(rgSet, type = "SnpI")$AddressB,
            getProbeInfo(rgSet, type = "SnpII")$AddressA)
    }
    indices <- which(rownames(rgSet) %in% addresses)

    # Filter RGChannelSet to retain only the selected probes
    rgSet[indices, ]
}

# Exported methods -------------------------------------------------------------

setMethod("show", signature(object = "RGChannelSet"), function(object) {
    callNextMethod()
    .show.annotation(annotation(object))
})

setMethod("annotation", signature(object = "RGChannelSet"), function(object) {
    object@annotation
})

setReplaceMethod(
    "annotation",
    signature(object = "RGChannelSet"),
    function(object, value) {
        object@annotation <- value
        object
    }
)

setMethod(
    "updateObject",
    signature(object = "RGChannelSet"),
    function(object, ..., verbose = FALSE) {
        if (verbose) message("updateObject(object = 'RGChannelSet')")
        if ("assayData" %in% names(getObjectSlots(object))) {
            # This is an ExpressionSet based object
            object <- RGChannelSet(
                Green = getObjectSlots(object)[["assayData"]][["Green"]],
                Red = getObjectSlots(object)[["assayData"]][["Red"]],
                colData = getObjectSlots(
                    getObjectSlots(object)[["phenoData"]])[["data"]],
                annotation = getObjectSlots(object)[["annotation"]])
        } else {
            object <- callNextMethod()
        }
        object
    }
)

setMethod("getManifest", signature(object = "RGChannelSet"), function(object) {
    maniString <- .getManifestString(object@annotation)
    if (!require(maniString, character.only = TRUE)) {
        stop(sprintf("cannot load manifest package %s", maniString))
    }
    updateObject(get(maniString))
})

setMethod("getBeta", signature(object = "RGChannelSet"), function(object, ...) {
    object <- preprocessRaw(object)
    callGeneric(object, ...)
})

setMethod(
    "mapToGenome",
    signature(object = "RGChannelSet"),
    function(object, ...) {
        object <- preprocessRaw(object)
        callGeneric(object, ...)
    }
)

setMethod(
    "combine",
    signature(x = "RGChannelSet", y = "RGChannelSet"),
    function(x, y, ...) {
        colDataFix <- .harmonizeDataFrames(
            .pDataFix(colData(x)),
            .pDataFix(colData(y)))
        colData(x) <- colDataFix$x
        colData(y) <- colDataFix$y
        cbind(x, y)
    }
)
