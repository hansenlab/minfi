setClass("RatioSet",
         representation(preprocessMethod = "character"),
         contains = "eSet",
         prototype = prototype(new("VersionedBiobase",
         versions = c(classVersion("eSet"), RatioSet = "1.0.0"))))

setValidity("RatioSet", function(object) {
    slotsPresent <- intersect(c("Beta", "M", "CN"), assayDataElementNames(object))
    msg <- NULL
    if(! any(c("Beta", "M") %in% slotsPresent))
        msg <- validMsg(msg, sprintf("objects of class '%s' needs to have assayData elements either 'Beta' or 'M' or both"),
                        class(object))
    msg <- validMsg(msg, assayDataValidMembers(assayData(object), slotsPresent))
    if (is.null(msg)) TRUE else msg
})

RatioSet <- function(Beta = NULL, M = NULL, CN=NULL, ...) {
    if(is.null(Beta) && is.null(M)) {
        Beta <- new("matrix")
    }
    if(!is.null(Beta) && is.null(M) && is.null(CN))
        Rset <- new("RatioSet", Beta = Beta, ...)
    if(!is.null(Beta) && is.null(M) && !is.null(CN))
        Rset <- new("RatioSet", Beta = Beta, CN = CN, ...)
    if(!is.null(Beta) && !is.null(M) && is.null(CN))
        Rset <- new("RatioSet", Beta = Beta, M = M, ...)
    if(!is.null(Beta) && !is.null(M) && !is.null(CN))
        Rset <- new("RatioSet", Beta = Beta, M = M, CN = CN, ...)
    if(is.null(Beta) && !is.null(M) && is.null(CN))
        Rset <- new("RatioSet", M = M, ...)
    if(is.null(Beta) && !is.null(M) && !is.null(CN))
        Rset <- new("RatioSet", M = M, CN = CN, ...)
    Rset
}

setMethod("show", "RatioSet", function(object) {
    .show.ExpressionSet(object)
    .show.annotation(annotation(object))
    .show.preprocessMethod(preprocessMethod(object))
})

setReplaceMethod("pData", signature(object = "RatioSet", value = "DataFrame"),
                 function(object, value) {
                     df <- as.data.frame(value)
                     pData(object) <- df
                     object
                 })

setMethod("getBeta", signature(object = "RatioSet"),
          function (object) {
              nms <- assayDataElementNames(object)
              if("Beta" %in% nms)
                  return(assayDataElement(object, "Beta"))
              if("M" %in% nms)
                  return(ilogit2(assayDataElement(object, "M")))
              stop("object does not contain either 'M' nor 'Beta' amongst assay slots")
          })

setMethod("getM", signature(object = "RatioSet"),
          function (object) {
              nms <- assayDataElementNames(object)
              if("M" %in% nms)
                  return(assayDataElement(object, "M"))
              if("Beta" %in% nms)
                  return(logit2(assayDataElement(object, "Beta")))
              stop("object does not contain either 'M' nor 'Beta' amongst assay slots")
          })

setMethod("getCN", signature(object = "RatioSet"),
          function (object) {
              nms <- assayDataElementNames(object)
              if("CN" %in% nms)
                  return(assayDataElement(object, "CN"))
              else
                  return(NULL)
          })

setMethod("preprocessMethod", signature(object = "RatioSet"),
          function(object) {
              object@preprocessMethod
          })

setMethod("mapToGenome", signature(object = "RatioSet"),
          function(object, drop = TRUE, mergeManifest = FALSE) {
              gr <- getLocations(object, mergeManifest = mergeManifest,
                                 orderByLocation = TRUE)
              object <- object[names(gr),]
              GenomicRatioSet(gr = gr, Beta = getBeta(object),
                              M = getM(object), CN = getCN(object),
                              pData = pData(object),
                              preprocessMethod = preprocessMethod(object),
                              annotation = annotation(object))
          })

setMethod("updateObject", signature(object = "RatioSet"),
          function(object, ..., verbose = FALSE) {
              if(object@annotation["annotation"] == "ilmn.v1.2")
                  object@annotation["annotation"] <- .default.450k.annotation
              object
          })
