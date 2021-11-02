# Exported classes -------------------------------------------------------------

setClass(
    "MethylSet",
    representation(preprocessMethod = "character", annotation = "character"),
    contains = "SummarizedExperiment"
)

# Validity methods -------------------------------------------------------------

setValidity("MethylSet", function(object) {
    msg <- validMsg(NULL, .checkAssayNames(object, c("Meth", "Unmeth")))
    if (is.null(msg)) TRUE else msg
})


# Exported functions -----------------------------------------------------------

MethylSet <- function(Meth = new("matrix"), Unmeth = new("matrix"),
                      annotation = "", preprocessMethod = "", ...) {
    # Check rownames, colnames
    assays <- SimpleList(Meth = Meth, Unmeth = Unmeth)
    new("MethylSet",
        SummarizedExperiment(assays = assays, ...),
        annotation = annotation,
        preprocessMethod = preprocessMethod
    )
}

dropMethylationLoci <- function(object, dropRS = TRUE, dropCH = TRUE) {
    stopifnot(class(object) %in% c(
        "MethylSet", "GenomicMethylSet", "RatioSet", "GenomicRatioSet"))
    dropRegEx <- ""
    if (dropRS) dropRegEx <- c(dropRegEx, "^rs")
    if (dropCH) dropRegEx <- c(dropRegEx, "^ch\\.")
    dropRegEx <- dropRegEx[nchar(dropRegEx) > 0]
    if (length(dropRegEx) == 0) return(object)
    dropRegEx <- sprintf("(%s)", paste(dropRegEx, collapse = "|"))
    whDrop <- grep(dropRegEx, rownames(object))
    if (length(whDrop) == 0) return(object)
    object[-whDrop, ]
}


# Exported methods -------------------------------------------------------------

setMethod("show", signature(object = "MethylSet"), function(object) {
    callNextMethod()
    .show.annotation(annotation(object))
    .show.preprocessMethod(preprocessMethod(object))
})

setMethod("getMeth", signature(object = "MethylSet"), function(object) {
    assay(object, "Meth")
})

setMethod("getUnmeth", signature(object = "MethylSet"), function(object) {
    assay(object, "Unmeth")
})

setMethod(
    "getBeta",
    signature(object = "MethylSet"),
    function(object, type = "", offset = 0, betaThreshold = 0) {
        if (type == "Illumina") offset <- 100
        .betaFromMethUnmeth(
            Meth = getMeth(object),
            Unmeth = getUnmeth(object),
            offset = offset,
            betaThreshold = betaThreshold)
    }
)

setMethod(
    "getM",
    signature(object = "MethylSet"),
    function(object, type = "", ...) {
        if (type == "") return(log2(getMeth(object) / getUnmeth(object)))
        if (type == "beta" || type == "Beta") logit2(getBeta(object, ...))
    }
)

setMethod("getCN", signature(object = "MethylSet"), function(object, ...) {
    log2(getMeth(object) + getUnmeth(object))
})

setMethod(
    "preprocessMethod",
    signature(object = "MethylSet"),
    function(object) object@preprocessMethod
)

setMethod("annotation", signature(object = "MethylSet"), function(object) {
    object@annotation
})

setReplaceMethod(
    "annotation",
    signature(object = "MethylSet"),
    function(object, value) {
        object@annotation <- value
        object
    }
)

setMethod(
    "getManifest",
    signature(object = "MethylSet"),
    function(object) {
        maniString <- .getManifestString(object@annotation)
        if (!require(maniString, character.only = TRUE))
            stop(sprintf("cannot load manifest package %s", maniString))
        get(maniString)
    }
)

setMethod(
    "mapToGenome",
    signature(object = "MethylSet"),
    function(object, mergeManifest = FALSE) {
        gr <- getLocations(
            object = object,
            mergeManifest = mergeManifest,
            orderByLocation = TRUE,
            lociNames = rownames(object))
        object <- object[names(gr), ]
        GenomicMethylSet(
            gr = gr,
            Meth = getMeth(object),
            Unmeth = getUnmeth(object),
            colData = colData(object),
            preprocessMethod = preprocessMethod(object),
            annotation = annotation(object),
            metadata = metadata(object))
    }
)

setMethod(
    "updateObject",
    signature(object = "MethylSet"),
    function(object, ..., verbose=FALSE) {
        if (verbose) message("updateObject(object = 'MethylSet')")
        if ("assayData" %in% names(getObjectSlots(object))) {
            # This is an ExpressionSet based object
            object <- MethylSet(
                Meth = getObjectSlots(object)[["assayData"]][["Meth"]],
                Unmeth = getObjectSlots(object)[["assayData"]][["Unmeth"]],
                colData = getObjectSlots(
                    getObjectSlots(object)[["phenoData"]])[["data"]],
                annotation = getObjectSlots(object)[["annotation"]],
                preprocessMethod = getObjectSlots(object)[["preprocessMethod"]])
        } else {
            object <- callNextMethod()
        }
        object
    }
)


setMethod(
    "ratioConvert",
    signature(object = "MethylSet"),
    function(object, what = c("beta", "M", "both"), keepCN = TRUE, ...) {
        what <- match.arg(what)
        if (what == "beta") {
            Beta <- getBeta(object, ...)
            M <- NULL
        }
        if (what == "M") {
            Beta <- NULL
            M <- getM(object, ...)
        }
        if (what == "both") {
            Beta <- getBeta(object, ...)
            M <- getM(object, ...)
        }
        CN <- if (keepCN) getCN(object) else NULL
        RatioSet(
            Beta = Beta,
            M = M,
            CN = CN,
            colData = colData(object),
            annotation = annotation(object),
            preprocessMethod = preprocessMethod(object),
            rowData = rowData(object),
            metadata = metadata(object))
    }
)

setMethod(
    "combine",
    signature(x = "MethylSet", y = "MethylSet"),
    function(x, y, ...) {
        colDataFix <- .harmonizeDataFrames(
            .pDataFix(colData(x)),
            .pDataFix(colData(y)))
        colData(x) <- colDataFix$x
        colData(y) <- colDataFix$y
        cbind(x, y)
    }
)
