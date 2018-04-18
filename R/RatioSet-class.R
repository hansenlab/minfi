# Exported classes -------------------------------------------------------------

setClass(
    "RatioSet",
    representation(preprocessMethod = "character", annotation = "character"),
    contains = "SummarizedExperiment")

# Validity methods -------------------------------------------------------------

setValidity("RatioSet", function(object) {
    msg <- validMsg(NULL, NULL)
    # NOTE: Do not check for the presence of 'CN' because it is optional
    msgBeta <- .checkAssayNames(object, c("Beta"))
    msgM <- .checkAssayNames(object, c("M"))
    if (!is.null(msgBeta) && !is.null(msgM)) {
        msg <- validMsg(
            msg,
            sprintf("objects of class '%s needs to have assays slots either 'Beta' or 'M' or both",
                    class(object)))
    }
    if (is.null(msg)) TRUE else msg
})

# Exported functions -----------------------------------------------------------

RatioSet <- function(Beta = NULL, M = NULL, CN = NULL, annotation = "",
                     preprocessMethod = "", ...) {
    msg <- "Need either 'Beta' or 'M' or both"
    if (is.null(Beta) && is.null(M)) stop(msg)
    assays <- SimpleList(Beta = Beta, M = M, CN = CN)
    assays <- assays[!vapply(assays, is.null, logical(1L))]
    new("RatioSet",
        SummarizedExperiment(assays = assays, ...),
        annotation = annotation,
        preprocessMethod = preprocessMethod
    )
}

# Exported methods -------------------------------------------------------------

setMethod("show", signature(object = "RatioSet"), function(object) {
    callNextMethod()
    .show.annotation(annotation(object))
    .show.preprocessMethod(preprocessMethod(object))
})

setMethod("preprocessMethod", signature(object = "RatioSet"), function(object) {
    object@preprocessMethod
})

setMethod("annotation", signature(object = "RatioSet"), function(object) {
    object@annotation
})

setReplaceMethod(
    "annotation",
    signature(object = "RatioSet"),
    function(object, value) {
        object@annotation <- value
        object
    }
)

setMethod("getBeta", signature(object = "RatioSet"), function(object) {
    nms <- names(assays(object, withDimnames = FALSE))
    if ("Beta" %in% nms) return(assay(object, "Beta"))
    if ("M" %in% nms) return(ilogit2(assay(object, "M")))
    stop("object does not contain either 'M' nor 'Beta' amongst assay slots")
})

setMethod("getM", signature(object = "RatioSet"), function(object) {
    nms <- names(assays(object, withDimnames = FALSE))
    if ("M" %in% nms) return(assay(object, "M"))
    if ("Beta" %in% nms) return(logit2(assay(object, "Beta")))
    stop("object does not contain either 'M' nor 'Beta' amongst assay slots")
})

setMethod("getCN", signature(object = "RatioSet"), function(object) {
    nms <- names(assays(object, withDimnames = FALSE))
    if ("CN" %in% nms) return(assay(object, "CN"))
    NULL
})

setMethod(
    "mapToGenome",
    signature(object = "RatioSet"),
    function(object, drop = TRUE, mergeManifest = FALSE) {
        gr <- getLocations(
            object = object,
            mergeManifest = mergeManifest,
            orderByLocation = TRUE)
        object <- object[names(gr),]
        GenomicRatioSet(
            gr = gr,
            Beta = getBeta(object),
            M = getM(object),
            CN = getCN(object),
            colData = colData(object),
            preprocessMethod = preprocessMethod(object),
            annotation = annotation(object))
    }
)

setMethod(
    "updateObject",
    signature(object = "RatioSet"),
    function(object, ..., verbose=FALSE) {
        if (verbose) message("updateObject(object = 'RatioSet')")
        if ("assayData" %in% names(getObjectSlots(object))) {
            # This is an ExpressionSet based object
            object <- RatioSet(
                Beta = getObjectSlots(object)[["assayData"]][["Beta"]],
                M = getObjectSlots(object)[["assayData"]][["M"]],
                CN = getObjectSlots(object)[["assayData"]][["CN"]],
                colData = getObjectSlots(
                    getObjectSlots(object)[["phenoData"]])[["data"]],
                annotation = getObjectSlots(object)[["annotation"]],
                preprocessMethod = getObjectSlots(object)[["preprocessMethod"]])
        }
        object
    }
)

setMethod(
    "combine",
    signature(x = "RatioSet", y = "RatioSet"),
    function(x, y, ...) {
        colDataFix <- .harmonizeDataFrames(.pDataFix(colData(x)),
                                           .pDataFix(colData(y)))
        colData(x) <- colDataFix$x
        colData(y) <- colDataFix$y
        cbind(x, y)
    }
)
