setMethod("sampleNames", signature(object = "RGChannelSet"),
          function(object) {
    colnames(object)
})
setMethod("sampleNames", signature(object = "MethylSet"),
          function(object) {
    colnames(object)
})
setMethod("sampleNames", signature(object = "RatioSet"),
          function(object) {
    colnames(object)
})
setMethod("sampleNames", signature(object = "GenomicMethylSet"),
          function(object) {
    colnames(object)
})
setMethod("sampleNames", signature(object = "GenomicRatioSet"),
          function(object) {
    colnames(object)
})

setReplaceMethod("sampleNames", signature(object = "RGChannelSet"),
          function(object, value) {
    colnames(object) <- value
    object
})
setReplaceMethod("sampleNames", signature(object = "MethylSet"),
          function(object, value) {
    colnames(object) <- value
    object
})
setReplaceMethod("sampleNames", signature(object = "RatioSet"),
          function(object, value) {
    colnames(object) <- value
    object
})
setReplaceMethod("sampleNames", signature(object = "GenomicMethylSet"),
          function(object, value) {
    colnames(object) <- value
    object
})
setReplaceMethod("sampleNames", signature(object = "GenomicRatioSet"),
          function(object, value) {
    colnames(object) <- value
    object
})

setMethod("featureNames", signature(object = "RGChannelSet"),
          function(object) {
    rownames(object)
})
setMethod("featureNames", signature(object = "MethylSet"),
          function(object) {
    rownames(object)
})
setMethod("featureNames", signature(object = "RatioSet"),
          function(object) {
    rownames(object)
})
setMethod("featureNames", signature(object = "GenomicMethylSet"),
          function(object) {
    rownames(object)
})
setMethod("featureNames", signature(object = "GenomicRatioSet"),
          function(object) {
    rownames(object)
})

setReplaceMethod("featureNames", signature(object = "RGChannelSet"),
          function(object, value) {
    rownames(object) <- value
    object
})
setReplaceMethod("featureNames", signature(object = "MethylSet"),
          function(object, value) {
    rownames(object) <- value
    object
})
setReplaceMethod("featureNames", signature(object = "RatioSet"),
          function(object, value) {
    rownames(object) <- value
    object
})
setReplaceMethod("featureNames", signature(object = "GenomicMethylSet"),
          function(object, value) {
    rownames(object) <- value
    object
})
setReplaceMethod("featureNames", signature(object = "GenomicRatioSet"),
          function(object, value) {
    rownames(object) <- value
    object
})

setMethod("pData", signature(object = "RGChannelSet"),
          function(object) {
    colData(object)
})
setMethod("pData", signature(object = "MethylSet"),
          function(object) {
    colData(object)
})
setMethod("pData", signature(object = "RatioSet"),
          function(object) {
    colData(object)
})
setMethod("pData", signature(object = "GenomicMethylSet"),
          function(object) {
    colData(object)
})
setMethod("pData", signature(object = "GenomicRatioSet"),
          function(object) {
    colData(object)
})

setReplaceMethod("pData", signature(object = "RGChannelSet", value = "DataFrame"),
                 function(object, value) {
    colData(object) <- value
    object
})
setReplaceMethod("pData", signature(object = "MethylSet", value = "DataFrame"),
                 function(object, value) {
    colData(object) <- value
    object
})
setReplaceMethod("pData", signature(object = "RatioSet", value = "DataFrame"),
                 function(object, value) {
    colData(object) <- value
    object
})
setReplaceMethod("pData", signature(object = "GenomicMethylSet", value = "DataFrame"),
                 function(object, value) {
    colData(object) <- value
    object
})
setReplaceMethod("pData", signature(object = "GenomicRatioSet", value = "DataFrame"),
                 function(object, value) {
    colData(object) <- value
    object
})

## setReplaceMethod("pData", signature(object = "GenomicMethylSet", value = "DataFrame"),
##                  function(object, value) {
##     object <- BiocGenerics:::replaceSlots(object, colData=value)
##     msg <- SummarizedExperiment:::.valid.SummarizedExperiment.assays_ncol(object)
##     if (!is.null(msg))
##         stop(msg)
##     object
## })
