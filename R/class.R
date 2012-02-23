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
              if (isCurrent(object)["RGChannelSet"]) return(object)
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
              if (isCurrent(object)["RGChannelSet"]) return(object)
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
    callNextMethod(object)
    preprocess <- object@preprocessMethod
    if(length(preprocess) == 0)
        preprocess <- c("unknown", "unknown", "unknown")
    cat("Preprocessing", preprocess[1],
        "minfi version", preprocess[2],
        "Manifest version", preprocess[3], fill = TRUE)
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

setClass("IlluminaMethylationManifest",
         representation(data = "environment",
                        annotation = "character"))

setValidity("IlluminaMethylationManifest", function(object) {
    msg <- NULL
    if(! "TypeI" %in% ls(object@data) || !is(object@data[["TypeI"]], "data.frame"))
        msg <- paste(msg, "'data' slot must contain a data.frame with TypeI probes", sep = "\n")
    if(! "TypeII" %in% ls(object@data) || !is(object@data[["TypeII"]], "data.frame"))
        msg <- paste(msg, "'data' slot must contain a data.frame with TypeII probes", sep = "\n")
    if(! "TypeControl" %in% ls(object@data) || !is(object@data[["TypeControl"]], "data.frame"))
        msg <- paste(msg, "'data' slot must contain a data.frame with Control probes", sep = "\n")
    if (is.null(msg)) TRUE else msg
})

setMethod("show", "IlluminaMethylationManifest", function(object) {
    cat("IlluminaMethylationManifest object\n")
    cat("Annotation:", object@annotation, "\n")
    cat("Number of type I probes:", nrow(object@data[["TypeI"]]), "\n")
    cat("Number of type II probes:", nrow(object@data[["TypeII"]]), "\n")
    cat("Number of control probes:", nrow(object@data[["TypeControl"]]), "\n")
})

manifestNew <- function(TypeI = new("data.frame"), TypeII = new("data.frame"),
                        TypeControl = new("data.frame"), annotation = "") {
    data <- new.env(parent = emptyenv())
    data[["TypeI"]] <- TypeI
    data[["TypeII"]] <- TypeII
    data[["TypeControl"]] <- TypeControl
    lockEnvironment(data, bindings = TRUE)
    manifest <- new("IlluminaMethylationManifest", annotation = annotation, data = data)
    manifest
}

