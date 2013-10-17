setClass("GenomicMethylSet",
         representation(preprocessMethod = "character", annotation = "character"),
         contains = "SummarizedExperiment")

setValidity("GenomicMethylSet", function(object) {
    msg <- validMsg(NULL, .checkAssayNames(object, c("Meth", "Unmeth")))
    if(class(rowData(object)) != "GRanges")
        msg <- validMsg(msg, sprintf("object of class '%s' needs to have a 'GRanges' in slot 'rowData'", class(object)))
    if (is.null(msg)) TRUE else msg
})

GenomicMethylSet <- function(gr = GRanges(), Meth = new("matrix"), Unmeth = new("matrix"),
                             pData = DataFrame(), annotation = "", preprocessMethod = "") {
    assays <- SimpleList(Meth = Meth, Unmeth = Unmeth)
    assays <- GenomicRanges:::.ShallowSimpleListAssays(data = assays)
    colData <- as(pData, "DataFrame")
    rowData <- as(gr, "GRanges")
    new("GenomicMethylSet", assays = assays, colData = colData,
        rowData = rowData, annotation = annotation, preprocessMethod = preprocessMethod)
}

setMethod("show", signature(object = "GenomicMethylSet"),
          function(object) {
              callNextMethod()
              .show.annotation(annotation(object))
              .show.preprocessMethod(preprocessMethod(object))
          })

setMethod("getMeth", signature(object = "GenomicMethylSet"),
          function(object) {
              assay(object, "Meth")
          })

setMethod("getUnmeth", signature(object = "GenomicMethylSet"),
          function(object) {
              assay(object, "Unmeth")
          })

setMethod("getBeta", signature(object = "GenomicMethylSet"),
          function(object, type = "", offset = 0, betaThreshold = 0) {
              if(type == "Illumina") {
                  offset <- 100
              }
              .betaFromMethUnmeth(Meth = getMeth(object), Unmeth = getUnmeth(object),
                                  offset = offset, betaThreshold = betaThreshold)
          })

setMethod("getM", signature(object = "GenomicMethylSet"),
          function (object, type = "", ...) {
              if(type == "")
                  return(log2(getMeth(object) / getUnmeth(object)))
              if(type == "beta" || type == "Beta")
                  return(logit2(getBeta(object, ...)))
          })

setMethod("getCN", signature(object = "GenomicMethylSet"),
          function(object, ...) {
              CN <- log2(getMeth(object) + getUnmeth(object))
              CN
          })

setMethod("mapToGenome", signature(object = "GenomicMethylSet"),
          function(object, ...) {
              object
          })

setMethod("preprocessMethod", signature(object = "GenomicMethylSet"),
          function(object) {
              object@preprocessMethod
          })

setMethod("annotation", signature(object = "GenomicMethylSet"),
          function(object) {
              object@annotation
          })

setMethod("pData", signature("GenomicMethylSet"),
          function(object) {
              colData(object)
          })

setReplaceMethod("pData", signature(object = "GenomicMethylSet", value = "DataFrame"),
                 function(object, value) {
                     object <- GenomicRanges:::clone(object, colData=value)
                     msg <- GenomicRanges:::.valid.SummarizedExperiment.colData_dims(object)
                     if (!is.null(msg))
                         stop(msg)
                     object
                 })

setMethod("sampleNames", signature("GenomicMethylSet"),
          function(object) {
              colnames(object)
          })

setMethod("featureNames", signature("GenomicMethylSet"),
          function(object) {
              rownames(object)
          })

setMethod("ratioConvert", signature(object = "GenomicMethylSet"),
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
              GenomicRatioSet(Beta = Beta, M = M, CN = CN,
                              gr = granges(object),
                              pData = pData(object),
                              annotation = annotation(object),
                              preprocessMethod = preprocessMethod(object))
          })

setMethod("updateObject", signature(object = "GenomicMethylSet"),
          function(object, ..., verbose = FALSE) {
              if(object@annotation["annotation"] == "ilmn.v1.2")
                  object@annotation["annotation"] <- .default.450k.annotation
              object
          })

