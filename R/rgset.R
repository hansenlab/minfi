setClass("RGChannelSet",
         contains = "eSet",
         prototype = prototype(new("VersionedBiobase",
         versions = c(classVersion("eSet"), RGChannelSet = "1.0.0"))))


setValidity("RGChannelSet", function(object) {
    msg <- validMsg(NULL, assayDataValidMembers(assayData(object),
                                                c("Red", "Green")))
    if (is.null(msg)) TRUE else msg
})

setMethod("initialize", "RGChannelSet",
          function(.Object, Red = new("matrix"), Green = new("matrix"), ...) {
              callNextMethod(.Object, Red = Red, Green = Green, ...)
})

RGChannelSet <- function(Green = new("matrix"), Red = new("matrix"), ...){
    rgSet <- new("RGChannelSet", Green = Green, Red = Red, ...)
    rgSet
}

setMethod("show", "RGChannelSet", function(object) {
    .show.ExpressionSet(object)
    .show.annotation(annotation(object))
})

setClass("RGChannelSetExtended",
         contains = "RGChannelSet")

setValidity("RGChannelSetExtended", function(object) {
    msg <- validMsg(NULL, assayDataValidMembers(assayData(object),
                                                c("Red", "Green", "RedSD", "GreenSD", "NBeads")))
    if (is.null(msg)) TRUE else msg
})

setMethod("initialize", "RGChannelSetExtended",
          function(.Object, Red = new("matrix"), Green = new("matrix"), RedSD = new("matrix"),
                   GreenSD = new("matrix"), NBeads = new("matrix"), ...) {
              callNextMethod(.Object, Red = Red, Green = Green, RedSD = RedSD,
                             GreenSD = GreenSD, NBeads = NBeads, ...)
})

RGChannelSetExtended <- function(Green = new("matrix"), Red = new("matrix"),
                                 GreenSD = new("matrix"), RedSD = new("matrix"),
                                 NBeads = new("matrix"), ...){
    rgSetEx <- new("RGChannelSet", Green = Green, Red = Red,
                   GreenSD = GreenSD, RedSD = RedSD, NBeads = NBeads, ...)
    rgSetEx
}

setMethod("updateObject", signature(object="RGChannelSet"),
          function(object, ..., verbose=FALSE) {
              if (verbose) message("updateObject(object = 'RGChannelSet')")
              ## object <- callNextMethod()
              if(object@annotation["annotation"] == "ilmn.v1.2")
                  object@annotation["annotation"] <- .default.450k.annotation
              if (isCurrent(object)["RGChannelSet"]) {
                  return(object)
              }
              if (! "RGChannelSet" %in% names(classVersion(object)))
                  newObject <- RGChannelSet(Red = assayDataElement(object, "Red"),
                                            Green = assayDataElement(object, "Green"),
                                            phenoData = updateObject(phenoData(object), ..., verbose=verbose),
                                            featureData = updateObject(featureData(object), ..., verbose=verbose),
                                            experimentData = updateObject(experimentData(object), ..., verbose=verbose),
                                            annotation = updateObject(annotation(object), ..., verbose=verbose),
                                            protocolData = updateObject(protocolData(object), ..., verbose=verbose))
              else
                  stop("unknown version update")
              newObject
          })

setMethod("updateObject", signature(object="RGChannelSetExtended"),
          function(object, ..., verbose=FALSE) {
              if (verbose) message("updateObject(object = 'RGChannelSetExtended')")
              ## object <- callNextMethod()
              if(object@annotation["annotation"] == "ilmn.v1.2")
                  object@annotation["annotation"] <- .default.450k.annotation
              if (isCurrent(object)["RGChannelSet"]) {
                  return(object)
              }
              if (! "RGChannelSet" %in% names(classVersion(object)))
                  newObject <- RGChannelSet(Red = assayDataElement(object, "Red"),
                                            Green = assayDataElement(object, "Green"),
                                            RedSD = assayDataElement(object, "RedSD"),
                                            GreenSD = assayDataElement(object, "GreenSD"),
                                            NBeads = assayDataElement(object, "NBeads"),
                                            phenoData = updateObject(phenoData(object), ..., verbose=verbose),
                                            featureData = updateObject(featureData(object), ..., verbose=verbose),
                                            experimentData = updateObject(experimentData(object), ..., verbose=verbose),
                                            annotation = updateObject(annotation(object), ..., verbose=verbose),
                                            protocolData = updateObject(protocolData(object), ..., verbose=verbose))
              else
                  stop("unknown version update")
              newObject
          })


getRed <- function(object) {
    .isRG(object)
    assayDataElement(object, "Red")
}

getGreen <- function(object) {
    .isRG(object)
    assayDataElement(object, "Green")
}

getOOB <- function(object) {
    .isRG(object)
    IRed   <- getProbeInfo(object, type = "I-Red")
    IGrn <- getProbeInfo(object, type = "I-Green")
    oob.green <- rbind(getGreen(object)[IRed$AddressA,], getGreen(object)[IRed$AddressB,])
    oob.red   <- rbind(getRed(object)[IGrn$AddressA,], getRed(object)[IGrn$AddressB,])
    return(list(Grn = oob.green, Red = oob.red))
}

getSnpBeta <- function(object){
    .isRG(object)
    manifest <- getManifest(object)
    snpProbesII <- getProbeInfo(manifest, type = "SnpII")$Name
    M.II <- getGreen(object)[getProbeInfo(object, type = "SnpII")$AddressA,]
    U.II <- getRed(object)[getProbeInfo(object, type = "SnpII")$AddressA,]
    snpProbesI.Green <- getProbeInfo(manifest, type = "SnpI")[getProbeInfo(manifest, type = "SnpI")$Color=="Grn",]
    snpProbesI.Red <- getProbeInfo(manifest, type = "SnpI")[getProbeInfo(manifest, type = "SnpI")$Color=="Red",]

    M.I.Red <- getRed(object)[snpProbesI.Red $AddressB,]
    U.I.Red <- getRed(object)[snpProbesI.Red $AddressA,]
    M.I.Green <- getGreen(object)[snpProbesI.Green$AddressB,]
    U.I.Green <- getGreen(object)[snpProbesI.Green$AddressA,]
    
    M <- rbind(M.II, M.I.Red, M.I.Green)
    U <- rbind(U.II, U.I.Red, U.I.Green)
    rownames(M) <- rownames(U) <- c(snpProbesII, snpProbesI.Red$Name, snpProbesI.Green$Name)
    beta <- M/(U+M+100)
    return(beta)                               
}

setMethod("getManifest", signature(object = "RGChannelSet"),
          function(object) {
              maniString <- .getManifestString(object@annotation)
              if(!require(maniString, character.only = TRUE))
                  stop(sprintf("cannot load manifest package %s", maniString))
              get(maniString)
          })

setMethod("getBeta", signature(object = "RGChannelSet"),
          function(object, ...) {
              object <- preprocessRaw(object)
              callGeneric(object)
          })

setMethod("mapToGenome", signature(object = "RGChannelSet"),
          function(object, ...) {
              object <- preprocessRaw(object)
              callGeneric(object)
          })
