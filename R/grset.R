setClass("GenomicRatioSet",
         representation(preprocessMethod = "character",
                        annotation = "character"),
         contains = "RangedSummarizedExperiment")

setValidity("GenomicRatioSet", function(object) {
    msg <- validMsg(NULL, NULL)
    ## It is intentional I do not check for the presence of 'CN'; it is optional
    msgBeta <- .checkAssayNames(object, c("Beta"))
    msgM <- .checkAssayNames(object, c("M"))
    if(!is.null(msgBeta) && !is.null(msgM))
        msg <- validMsg(msg, sprintf("objects of class '%s needs to have assays slots either 'Beta' or 'M' or both",
                                     class(object)))
    if(!is(rowRanges(object), "GRanges"))
        msg <- validMsg(msg, sprintf("object of class '%s' needs to have a 'GRanges' in slot 'rowRanges'", class(object)))
    if (is.null(msg)) TRUE else msg
})


GenomicRatioSet <- function(gr = GRanges(), Beta = NULL, M = NULL, CN = NULL,
                            annotation = "", preprocessMethod = "", ...) {
    msg <- "Need either 'Beta' or 'M' or both"
    if(is.null(Beta) && is.null(M))
        stop(msg)
    assays <- SimpleList(Beta = Beta, M = M, CN = CN)
    assays <- assays[!sapply(assays, is.null)]
    new("GenomicRatioSet",
        SummarizedExperiment(
            assays = assays,
            rowRanges = as(gr, "GRanges"),
            ...),
        annotation = annotation,
        preprocessMethod = preprocessMethod
    )
}

setMethod("show", signature(object = "GenomicRatioSet"),
          function(object) {
              callNextMethod()
              .show.annotation(annotation(object))
              .show.preprocessMethod(preprocessMethod(object))
          })

setMethod("getBeta", signature(object = "GenomicRatioSet"),
          function (object) {
              nms <- names(assays(object, withDimnames = FALSE))
              if("Beta" %in% nms)
                  return(assay(object, "Beta"))
              if("M" %in% nms)
                  return(ilogit2(assay(object, "M")))
              stop("object does not contain either 'M' nor 'Beta' amongst assay slots")
          })

setMethod("getM", signature(object = "GenomicRatioSet"),
          function (object) {
              nms <- names(assays(object, withDimnames = FALSE))
              if("M" %in% nms)
                  return(assay(object, "M"))
              if("Beta" %in% nms)
                  return(logit2(assay(object, "Beta")))
              stop("object does not contain either 'M' nor 'Beta' amongst assay slots")

          })

setMethod("getCN", signature(object = "GenomicRatioSet"),
          function (object) {
              nms <- names(assays(object, withDimnames = FALSE))
              if("CN" %in% nms)
                  return(assay(object, "CN"))
              return(NULL)
          })

# FIXME: enforce something a bit less lackadaisical?
setMethod("getSNPs", signature(object = "GenomicRatioSet"),
          function (object) {
              if (!all(colnames(object) %in% colnames(metadata(object)$SNPs))) {
                missed <- setdiff(colnames(object), 
                                  colnames(metadata(object)$SNPs))
                stop("Error: columns for ", paste(missed, sep=", "),
                     " are missing from metadata(object)$SNPs.")
              }
              if ("SNPs" %in% names(metadata(object))) { 
                return(metadata(object)$SNPs[, colnames(object)])
              } else {
                message("No SNPs found in your GenomicRatioSet's metadata.")
                return(NULL)
              }
          })

setMethod("mapToGenome", signature(object = "GenomicRatioSet"),
          function(object, ...) {
              object
          })

setMethod("preprocessMethod", signature(object = "GenomicRatioSet"),
          function(object) {
              object@preprocessMethod
          })

setMethod("annotation", signature(object = "GenomicRatioSet"),
          function(object) {
              object@annotation
          })

setReplaceMethod("annotation", signature(object = "GenomicRatioSet"),
          function(object, value) {
    object@annotation <- value
    object
})

setMethod("updateObject", signature(object = "GenomicRatioSet"),
          function(object, ..., verbose = FALSE) {
              if(object@annotation["annotation"] == "ilmn.v1.2")
                  object@annotation["annotation"] <- .default.450k.annotation
              object
          })

setMethod("combine", signature(x = "GenomicRatioSet", y = "GenomicRatioSet"),
          function(x, y, ...) {
    colDataFix <- .harmonizeDataFrames(.pDataFix(colData(x)),
                                       .pDataFix(colData(y)))
    colData(x) <- colDataFix$x
    colData(y) <- colDataFix$y
    cbind(x,y)
})
