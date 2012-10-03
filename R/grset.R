setClass("GenomicRatioSet",
         representation(preprocessMethod = "character",
                        annotation = "character"),
         contains = "SummarizedExperiment")

setValidity("GenomicRatioSet", function(object) {
    msg <- validMsg(NULL, .checkAssayNames(object, c("CN")))
    msgBeta <- .checkAssayNames(object, c("Beta"))
    msgM <- .checkAssayNames(object, c("M"))
    if(!is.null(msgBeta) && !is.null(msgM))
        msg <- validMsg(msg, sprintf("objects of class '%s% needs to have assays slots either 'Beta' or 'M' or both",
                                     class(object)))
    if(class(rowData(object)) != "GRanges")
        msg <- validMsg(msg, sprintf("object of class '%s' needs to have a 'GRanges' in slot 'rowData'", class(object)))
    if (is.null(msg)) TRUE else msg
})


GenomicRatioSet <- function(gr, Beta = NULL, M = NULL, CN = NULL, pData, annotation, preprocessMethod) {
    msg <- "Need either 'Beta' or 'M' or both"
    if(is.null(Beta) && is.null(M))
        stop(msg)
    if(!is.null(Beta) && is.null(M))
        assays <- SimpleList(Beta = Beta, CN = CN)
    if(is.null(Beta) && !is.null(M))
        assays <- SimpleList(M = M, CN = CN)
    if(!is.null(Beta) && !is.null(M))
        assays <- SimpleList(M = M, Beta = Beta, CN = CN)
    assays <- GenomicRanges:::.ShallowSimpleListAssays$new(data = assays)
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
                  return(logit2(assay(object, "M")))
              stop("object does not contain either 'M' nor 'Beta' amongst assay slots")

          })



setGeneric("getCN", function(object, ...) standardGeneric("getCN"))

setMethod("getCN", signature(object = "GenomicRatioSet"),
          function (object) {
            nms <- names(assays(object, withDimnames = FALSE))
            if("CN" %in% nms)
                return(assay(object, "CN"))
            stop("object does not contain either 'CN' amongst assay slots")
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
