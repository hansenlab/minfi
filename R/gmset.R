setClass("GenomicMethylSet",
         representation(preprocessMethod = "character", annotation = "character"),
         contains = "SummarizedExperiment")

## FIXME: add show method

setValidity("GenomicMethylSet", function(object) {
    msg <- validMsg(NULL, .checkAssayNames(object, c("Meth", "Unmeth")))
    if(class(rowData(object)) != "GRanges")
        msg <- validMsg(msg, sprintf("object of class '%s' needs to have a 'GRanges' in slot 'rowData'", class(object)))
    if (is.null(msg)) TRUE else msg
})

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
          function (object, type = "", ...) {
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

setReplaceMethod("pData", c("GenomicMethylSet", "DataFrame"),
    function(x, ..., value)
{
    x <- clone(x, ..., colData=value)
    msg <- GenomicRanges:::.valid.SummarizedExperiment.colData_dims(x)
    if (!is.null(msg))
        stop(msg)
    x
})

setMethod("sampleNames", signature("GenomicMethylSet"),
          function(object) {
              colnames(object)
          })

setMethod("featureNames", signature("GenomicMethylSet"),
          function(object) {
              rownames(object)
          })

