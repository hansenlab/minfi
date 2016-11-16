setClass("MethylSet",
         representation(preprocessMethod = "character", annotation = "character"),
         contains = "SummarizedExperiment")

setValidity("MethylSet", function(object) {
    msg <- validMsg(NULL, .checkAssayNames(object, c("Meth", "Unmeth")))
    if(is.null(msg)) TRUE else msg
})

MethylSet <- function(Meth = new("matrix"), Unmeth = new("matrix"),
                      pData = DataFrame(), annotation = "", preprocessMethod = "",
                      rowData = NULL, metadata = list()) {
    ## Check rownames, colnames
    assays <- SimpleList(Meth = Meth, Unmeth = Unmeth)
    new("MethylSet",
        SummarizedExperiment(
            assays = assays,
            colData = as(pData, "DataFrame"),
            rowData = rowData,
            metadata = metadata
        ),
        annotation = annotation,
        preprocessMethod = preprocessMethod
        )
}

setMethod("show", signature(object = "MethylSet"),
          function(object) {
    callNextMethod()
    .show.annotation(annotation(object))
    .show.preprocessMethod(preprocessMethod(object))
})

setMethod("pData", signature("MethylSet"),
          function(object) {
    colData(object)
})

setReplaceMethod("pData", signature(object = "MethylSet", value = "DataFrame"),
                 function(object, value) {
    ## FIXME: Different from the GenomicMethylSet function; which one is right?
    ## This one is easier
    colData(object) <- value
    object
})

setMethod("getMeth", signature(object = "MethylSet"),
          function(object) {
    assay(object, "Meth")
})

setMethod("getUnmeth", signature(object = "MethylSet"),
          function(object) {
    assay(object, "Unmeth")
})

setMethod("getBeta", signature(object = "MethylSet"),
          function(object, type = "", offset = 0, betaThreshold = 0) {
    if(type == "Illumina") {
        offset <- 100
    }
    .betaFromMethUnmeth(Meth = getMeth(object), Unmeth = getUnmeth(object),
                        offset = offset, betaThreshold = betaThreshold)
})

setMethod("getM", signature(object = "MethylSet"),
          function (object, type = "", ...) {
    if(type == "")
        return(log2(getMeth(object) / getUnmeth(object)))
    if(type == "beta" || type == "Beta")
        return(logit2(getBeta(object, ...)))
})

setMethod("getCN", signature(object = "MethylSet"),
          function(object, ...) {
    CN <- log2(getMeth(object) + getUnmeth(object))
    CN
})

setMethod("preprocessMethod", signature(object = "MethylSet"),
          function(object) {
    object@preprocessMethod
})

setMethod("annotation", signature(object = "MethylSet"),
          function(object) {
    object@annotation
})

setMethod("sampleNames", signature("MethylSet"),
          function(object) {
    colnames(object)
})

setMethod("featureNames", signature("MethylSet"),
          function(object) {
    rownames(object)
})

setMethod("getManifest", signature(object = "MethylSet"),
          function(object) {
    maniString <- .getManifestString(object@annotation)
    if(!require(maniString, character.only = TRUE))
        stop(sprintf("cannot load manifest package %s", maniString))
    get(maniString)
})

setMethod("mapToGenome", signature(object = "MethylSet"),
          function(object, mergeManifest = FALSE) {
    gr <- getLocations(object, mergeManifest = mergeManifest,
                       orderByLocation = TRUE, lociNames = featureNames(object))
    object <- object[names(gr),]
    GenomicMethylSet(gr = gr, Meth = getMeth(object),
                     Unmeth = getUnmeth(object),
                     pData = pData(object),
                     preprocessMethod = preprocessMethod(object),
                     annotation = annotation(object),
                     metadata = metadata(object))
})

setMethod("updateObject", signature(object = "MethylSet"),
          function(object, ..., verbose=FALSE) {
    if (verbose) message("updateObject(object = 'MethylSet')")
    if("assayData" %in% names(getObjectSlots(object))) {
        ## This is an ExpressionSet based object
        newObject <- MethylSet(Meth = getObjectSlots(object)[["assayData"]][["Meth"]],
                               Unmeth = getObjectSlots(object)[["assayData"]][["Unmeth"]],
                               pData = getObjectSlots(getObjectSlots(object)[["phenoData"]])[["data"]],
                               annotation = getObjectSlots(object)[["annotation"]],
                               preprocessMethod = getObjectSlots(object)[["preprocessMethod"]])
    }
    newObject
})


setMethod("ratioConvert", signature(object = "MethylSet"),
          function(object, what = c("beta", "M", "both"), keepCN = TRUE, ...) {
    what <- match.arg(what)
    if(what == "beta") {
        Beta <- getBeta(object, ...)
        M <- NULL
    }
    if(what == "M") {
        Beta <- NULL
        M <- getM(object, ...)
    }
    if(what == "both") {
        Beta <- getBeta(object, ...)
        M <- getM(object, ...)
    }
    if(keepCN) {
        CN <- getCN(object)
    } else {
        CN <- NULL
    }
    RatioSet(Beta = Beta, M = M, CN = CN,
             pData = colData(object),
             annotation = annotation(object),
             preprocessMethod = preprocessMethod(object),
             rowData = rowData(object),
             metadata = metadata(object))
})

dropMethylationLoci <- function(object, dropRS = TRUE, dropCH = TRUE) {
    stopifnot(class(object) %in% c("MethylSet", "GenomicMethylSet", "RatioSet", "GenomicRatioSet"))
    dropRegEx <- ""
    if(dropRS)
        dropRegEx <- c(dropRegEx, "^rs")
    if(dropCH)
        dropRegEx <- c(dropRegEx, "^ch\\.")
    dropRegEx <- dropRegEx[nchar(dropRegEx) > 0]
    if(length(dropRegEx) == 0) return(object)
    dropRegEx <- sprintf("(%s)", paste(dropRegEx, collapse = "|"))
    whDrop <- grep(dropRegEx, featureNames(object))
    if(length(whDrop) == 0)
        return(object)
    object[-whDrop, ]
}

setMethod("combine", signature(x = "MethylSet", y = "MethylSet"),
          function(x, y, ...) {
    pData(x) <- .pDataFix(pData(x))
    pData(y) <- .pDataFix(pData(y))
    callNextMethod()
})
