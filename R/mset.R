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


setMethod("getMeth", signature(object = "MethylSet"),
          function(object) {
              assayDataElement(object, "Meth")
          })

setMethod("getUnmeth", signature(object = "MethylSet"),
          function(object) {
              assayDataElement(object, "Unmeth")
          })

setMethod("getBeta", signature(object = "MethylSet"),
          function(object, type = c("minfi", "Illumina")) {
              type <- match.arg(type)
              Meth <- pmax(assayDataElement(object, "Meth"), 0)
              Unmeth <- pmax(assayDataElement(object, "Unmeth"), 0)
              if(type == "minfi")
                  return(Meth / (Meth + Unmeth))
              if(type == "Illumina")
                  return(Meth / (Meth + Unmeth + 100))
          })

setMethod("getM", signature(object = "MethylSet"),
          function (object, type = c("minfi", "Illumina")) {  
              type <- match.arg(type)
              beta <- getBeta(object, type = type)
              beta <- pmin(beta, 0.999)
              beta <- pmax(beta, 0.001)
              logit2(beta)
          })
