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

MethylSet <- function(Meth, Unmeth, phenoData, annotation = NULL) {
    if(is.null(annotation))
        annotation <- c(array = "IlluminaHumanMethylation450k",
                        annotation = .default.450k.annotation)
    stopifnot(rownames(Meth) == rownames(Unmeth))
    stopifnot(colnames(Meth) == colnames(Unmeth))
    stopifnot(colnames(Meth) == rownames(phenoData))
    tmp <- matrix(nrow = 0, ncol = ncol(Meth))
    out <- new("MethylSet", Meth = tmp, Unmeth = tmp,
               phenoData = phenoData, annotation = annotation)
    assayDataElement(out, "Meth") <- Meth
    assayDataElement(out, "Unmeth") <- Unmeth
    featureData(out) <- AnnotatedDataFrame(data = data.frame(row.names = row.names(Meth)),
                                           dimLabels = c("featureNames", "featureColumns"))
    out
}

setMethod("show", "MethylSet", function(object) {
    .show.ExpressionSet(object)
    .show.annotation(annotation(object))
    .show.preprocessMethod(preprocessMethod(object))
})

setReplaceMethod("pData", signature(object = "MethylSet", value = "DataFrame"),
                 function(object, value) {
                     df <- as.data.frame(value)
                     pData(object) <- df
                     object
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

setMethod("getCN", signature(object = "MethylSet"),
          function(object, ...) {
              CN <- log2(getMeth(object) + getUnmeth(object))
              CN
          })

setMethod("preprocessMethod", signature(object = "MethylSet"),
          function(object) {
              object@preprocessMethod
          })

setMethod("getManifest", signature(object = "MethylSet"),
          function(object) {
              maniString <- .getManifestString(object@annotation)
              if(!require(maniString, character.only = TRUE))
                  stop(sprintf("cannot load manifest package %s", maniString))
              get(maniString)
          })

setMethod("mapToGenome", signature(object = "MethylSet"),
          function(object, mergeManifest = FALSE) {
              gr <- getLocations(object, mergeManifest = mergeManifest,
                                 orderByLocation = TRUE, lociNames = featureNames(object))
              object <- object[names(gr),]
              GenomicMethylSet(gr = gr, Meth = getMeth(object),
                               Unmeth = getUnmeth(object),
                               pData = pData(object),
                               preprocessMethod = preprocessMethod(object),
                               annotation = annotation(object))
          })

setMethod("updateObject", signature(object = "MethylSet"),
          function(object, ..., verbose=FALSE) {
              if (verbose) message("updateObject(object = 'MethylSet')")
              ## object <- callNextMethod()
              if (isCurrent(object)["MethylSet"]) {
                  if(object@annotation["annotation"] == "ilmn.v1.2")
                      object@annotation["annotation"] <- .default.450k.annotation
                  return(object)
              }
              if (! "MethylSet" %in% names(classVersion(object)))
                  newObject <- new("MethylSet", preprocessMethod = object@preprocessMethod,
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

setMethod("ratioConvert", signature(object = "MethylSet"),
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
              RatioSet(Beta = Beta, M = M, CN = CN,
                       phenoData = phenoData(object),
                       annotation = annotation(object),
                       featureData = featureData(object),
                       preprocessMethod = preprocessMethod(object))
          })
    


dropMethylationLoci <- function(object, dropRS = TRUE, dropCH = TRUE) {
    stopifnot(class(object) %in% c("MethylSet", "GenomicMethylSet", "RatioSet", "GenomicRatioSet"))
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
