setClass("MethylSet",
         representation(preprocessMethod = "character"),
         contains = "eSet",
         prototype = prototype(new("VersionedBiobase",
         versions = c(classVersion("eSet"), MethylSet = "1.0.0"))))

setValidity("MethylSet", function(object) {
    msg <- validMsg(NULL, assayDataValidMembers(assayData(object),
                                                c("Meth", "Unmeth")))
    if (is.null(msg)) TRUE else msg
})

MethylSet <- function(Meth = new("matrix"), Unmeth = new("matrix"),  ...) {
    Mset <- new("MethylSet", Meth = Meth, Unmeth = Unmeth, ...)
    Mset
}

setMethod("show", "MethylSet", function(object) {
    .show.ExpressionSet(object)
    .show.annotation(object@annotation)
    .show.preprocessMethod(object@preprocessMethod)
})

setMethod("getMeth", signature(object = "MethylSet"),
          function(object) {
              assayDataElement(object, "Meth")
          })

setMethod("getUnmeth", signature(object = "MethylSet"),
          function(object) {
              assayDataElement(object, "Unmeth")
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

setMethod("getManifest", signature(object = "MethylSet"),
          function(object) {
              maniString <- .getManifestString(object@annotation)
              if(!require(maniString, character.only = TRUE))
                  stop(sprintf("cannot load manifest package %s", maniString))
              get(maniString)
          })

setMethod("updateObject", signature(object="MethylSet"),
          function(object, ..., verbose=FALSE) {
              if (verbose) message("updateObject(object = 'MethylSet')")
              ## object <- callNextMethod()
              if (isCurrent(object)["MethylSet"]) return(object)
              if (! "MethylSet" %in% names(classVersion(object)))
                  newObject <- MethylSet(preprocessMethod = object@preprocessMethod,
                                         Meth = assayDataElement(object, "Meth"),
                                         Unmeth = assayDataElement(object, "Unmeth"),
                                         phenoData = updateObject(phenoData(object), ..., verbose=verbose),
                                         featureData = updateObject(featureData(object), ..., verbose=verbose),
                                         experimentData = updateObject(experimentData(object), ..., verbose=verbose),
                                         annotation = updateObject(annotation(object), ..., verbose=verbose),
                                         protocolData = updateObject(protocolData(object), ..., verbose=verbose))
              else
                  stop("unknown version update")
              newObject
          })

dropMethylationLoci <- function(object, dropRS = TRUE, dropCH = TRUE) {
    stopifnot(class(object) %in% c("MethylSet", "GenomicMethylSet"))
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
