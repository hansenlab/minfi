setClass("GenomicMethylSet",
         representation(preprocessMethod = "character", annotation = "character"),
         contains = "SummarizedExperiment")

## FIXME: add show method
## FIXME: add validity, check for GRanges, not GRangesList, meth, unmeth

GenomicMethylSet <- function(gr, Meth, Unmeth, pData, annotation, preprocessMethod) {
    assays <- SimpleList(Meth = Meth, Unmeth = Unmeth)
    assays <- GenomicRanges:::.ShallowSimpleListAssays$new(data = assays)
    colData <- as(pData, "DataFrame")
    rowData <- as(gr, "GRanges")
    new("GenomicMethylSet", assays = assays, colData = colData,
        rowData = rowData, annotation = annotation, preprocessMethod = preprocessMethod)
}

setMethod("show", signature(object = "GenomicMethylSet"),
          function(object) {
              callNextMethod()
              minfi:::.show.annotation(annotation(object))
              minfi:::.show.preprocessMethod(preprocessMethod(object))
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
          function (object, ...) {
              if(type == "")
                  return(log2(getMeth(object) / getUnmeth(object)))
              if(type == "beta" || type == "Beta")
                  return(logit2(getBeta(object, ...)))
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

setMethod("sampleNames", signature("GenomicMethylSet"),
          function(object) {
              colnames(object)
          })

setMethod("featureNames", signature("GenomicMethylSet"),
          function(object) {
              rownames(object)
          })

