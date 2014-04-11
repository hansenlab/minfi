setClass("GenomicRatioSet",
         representation(preprocessMethod = "character",
                        annotation = "character"),
         contains = "SummarizedExperiment")

setValidity("GenomicRatioSet", function(object) {
    msg <- validMsg(NULL, NULL)
    ## It is intentional I do not check for the presence of 'CN'; it is optional
    msgBeta <- .checkAssayNames(object, c("Beta"))
    msgM <- .checkAssayNames(object, c("M"))
    if(!is.null(msgBeta) && !is.null(msgM))
        msg <- validMsg(msg, sprintf("objects of class '%s needs to have assays slots either 'Beta' or 'M' or both",
                                     class(object)))
    if(!is(rowData(object), "GRanges"))
        msg <- validMsg(msg, sprintf("object of class '%s' needs to have a 'GRanges' in slot 'rowData'", class(object)))
    if (is.null(msg)) TRUE else msg
})


GenomicRatioSet <- function(gr = GRanges(), Beta = NULL, M = NULL, CN = NULL,
                            pData = DataFrame(),
                            annotation = "",
                            preprocessMethod = "") {
    msg <- "Need either 'Beta' or 'M' or both"
    if(is.null(Beta) && is.null(M))
        stop(msg)
    assays <- SimpleList(Beta = Beta, M = M, CN = CN)
    assays <- GenomicRanges:::.ShallowSimpleListAssays(data = assays[!sapply(assays, is.null)])
    colData <- as(pData, "DataFrame")
    rowData <- as(gr, "GRanges")
    new("GenomicRatioSet", assays = assays, colData = colData,
        rowData = rowData, annotation = annotation, preprocessMethod = preprocessMethod)
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

setMethod("pData", signature("GenomicRatioSet"),
          function(object) {
              colData(object)
          })

setReplaceMethod("pData", c("GenomicRatioSet", "DataFrame"), function(object, value) {
    object <- GenomicRanges:::clone(object, colData=value)
    msg <- GenomicRanges:::.valid.SummarizedExperiment.colData_dims(object)
    if (!is.null(msg))
        stop(msg)
    object
})

setMethod("mapToGenome", signature(object = "GenomicRatioSet"),
          function(object, ...) {
              object
          })

setMethod("sampleNames", signature("GenomicRatioSet"),
          function(object) {
              colnames(object)
          })

setMethod("featureNames", signature("GenomicRatioSet"),
          function(object) {
              rownames(object)
          })

setMethod("preprocessMethod", signature(object = "GenomicRatioSet"),
          function(object) {
              object@preprocessMethod
          })

setMethod("annotation", signature(object = "GenomicRatioSet"),
          function(object) {
              object@annotation
          })

setMethod("updateObject", signature(object = "GenomicRatioSet"),
          function(object, ..., verbose = FALSE) {
              if(object@annotation["annotation"] == "ilmn.v1.2")
                  object@annotation["annotation"] <- .default.450k.annotation
              object
          })
