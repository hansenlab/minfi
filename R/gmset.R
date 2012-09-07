setClass("GenomicMethylSet",
         representation(preprocessMethod = "character"),
         contains = "SummarizedExperiment")

## FIXME: add show method
## FIXME: add validity, check for GRanges, not GRangesList, meth, unmeth

GenomicMethylSet <- function(gr, Meth, Unmeth, pData) {
    assays <- SimpleList(Meth = Meth, Unmeth = Unmeth)
    assays <- GenomicRanges:::.ShallowSimpleListAssays$new(data = assays)
    colData <- as(pData, "DataFrame")
    rowData <- as(gr, "GRanges")
    new("GenomicMethylSet", assays = assays, colData = colData,
        rowData = rowData)
}

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

setMethod("granges", signature("GenomicMethylSet"),
          function(x, ...) {
              rowData(x)
          })

setMethod("start", signature("GenomicMethylSet"),
          function(x, ...) {
              start(rowData(x))
          })

setMethod("end", signature("GenomicMethylSet"),
          function(x, ...) {
              end(rowData(x))
          })

setMethod("width", signature("GenomicMethylSet"),
          function(x) {
              width(rowData(x))
          })

setMethod("strand", signature("GenomicMethylSet"),
          function(x, ...) {
              strand(rowData(x))
          })

setMethod("seqnames", signature("GenomicMethylSet"),
          function(x) {
              seqnames(rowData(x))
          })

setMethod("seqlevels", signature("GenomicMethylSet"),
          function(x) {
              seqlevels(rowData(x))
          })

setMethod("seqlengths", signature("GenomicMethylSet"),
          function(x) {
              seqlengths(rowData(x))
          })

setMethod("genome", signature("GenomicMethylSet"),
          function(x) {
              genome(rowData(x))
          })


