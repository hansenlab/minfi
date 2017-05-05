setClass("RGChannelSet",
         representation(annotation = "character"),
         contains = "SummarizedExperiment")

setValidity("RGChannelSet", function(object) {
    msg <- validMsg(NULL, .checkAssayNames(object, c("Red", "Green")))
    if(is.null(msg)) TRUE else msg
})

RGChannelSet <- function(Green = new("matrix"), Red = new("matrix"),
                         annotation = "", ...) {
    ## Check rownames, colnames
    assays <- SimpleList(Green = Green, Red = Red)
    new("RGChannelSet",
        SummarizedExperiment(assays = assays, ...),
        annotation = annotation
        )
}

setMethod("show", signature(object = "RGChannelSet"),
          function(object) {
    callNextMethod()
    .show.annotation(annotation(object))
})

setMethod("annotation", signature(object = "RGChannelSet"),
          function(object) {
    object@annotation
})

setReplaceMethod("annotation", signature(object = "RGChannelSet"),
          function(object, value) {
    object@annotation <- value
    object
})

setClass("RGChannelSetExtended",
         contains = "RGChannelSet")

setValidity("RGChannelSetExtended", function(object) {
    msg <- validMsg(NULL, .checkAssayNames(object,
                                           c("Red", "Green", "RedSD", "GreenSD", "NBeads")))
    if (is.null(msg)) TRUE else msg
})

RGChannelSetExtended <- function(Green = new("matrix"), Red = new("matrix"),
                                 GreenSD = new("matrix"), RedSD = new("matrix"),
                                 NBeads = new("matrix"), annotation = "",
                                 ...) {
    ## Check rownames, colnames
    assays <- SimpleList(Green = Green, Red = Red,
                         GreenSD = GreenSD, RedSD = RedSD,
                         NBeads = NBeads)
    new("RGChannelSetExtended",
        SummarizedExperiment(assays = assays, ...),
        annotation = annotation
        )
}

setMethod("updateObject", signature(object = "RGChannelSet"),
          function(object, ..., verbose=FALSE) {
    if (verbose) message("updateObject(object = 'RGChannelSet')")
    if("assayData" %in% names(getObjectSlots(object))) {
        ## This is an ExpressionSet based object
        object <- RGChannelSet(Green = getObjectSlots(object)[["assayData"]][["Green"]],
                               Red = getObjectSlots(object)[["assayData"]][["Red"]],
                               colData = getObjectSlots(getObjectSlots(object)[["phenoData"]])[["data"]],
                               annotation = getObjectSlots(object)[["annotation"]])
    }
    object
})

setMethod("updateObject", signature(object = "RGChannelSetExtended"),
          function(object, ..., verbose=FALSE) {
    if (verbose) message("updateObject(object = 'RGChannelSetExtended')")
    if("assayData" %in% names(getObjectSlots(object))) {
        ## This is an ExpressionSet based object
        object <- RGChannelSetExtended(Green = getObjectSlots(object)[["assayData"]][["Green"]],
                                       Red = getObjectSlots(object)[["assayData"]][["Red"]],
                                       GreenSD = getObjectSlots(object)[["assayData"]][["GreenSD"]],
                                       RedSD = getObjectSlots(object)[["assayData"]][["RedSD"]],
                                       NBeads = getObjectSlots(object)[["assayData"]][["NBeads"]],
                                       colData = getObjectSlots(getObjectSlots(object)[["phenoData"]])[["data"]],
                                       annotation = getObjectSlots(object)[["annotation"]])
    }
    object
})


getRed <- function(object) {
    .isRGOrStop(object)
    assay(object, "Red")
}

getGreen <- function(object) {
    .isRGOrStop(object)
    assay(object, "Green")
}

getOOB <- function(object) {
    .isRGOrStop(object)
    IRed   <- getProbeInfo(object, type = "I-Red")
    IGrn <- getProbeInfo(object, type = "I-Green")
    oob.green <- rbind(getGreen(object)[IRed$AddressA,,drop=FALSE], getGreen(object)[IRed$AddressB,,drop=FALSE])
    oob.red   <- rbind(getRed(object)[IGrn$AddressA,,drop=FALSE], getRed(object)[IGrn$AddressB,,drop=FALSE])
    return(list(Grn = oob.green, Red = oob.red))
}

getNBeads <- function(object) {
    if(!is(object, "RGChannelSetExtended"))
        stop(sprintf("object is of class '%s', but needs to be of class 'RGChannelSetExtended'", class(object)))
    assay(object, "NBeads")
}


getSnpBeta <- function(object){
    .isRGOrStop(object)

    snpProbesII <- getProbeInfo(object, type = "SnpII")$Name
    M.II <- getGreen(object)[getProbeInfo(object, type = "SnpII")$AddressA,,drop=FALSE]
    U.II <- getRed(object)[getProbeInfo(object, type = "SnpII")$AddressA,,drop=FALSE]

    snpProbesI.Green <- getProbeInfo(object, type = "SnpI")[getProbeInfo(object, type = "SnpI")$Color=="Grn",,drop=FALSE]
    snpProbesI.Red <- getProbeInfo(object, type = "SnpI")[getProbeInfo(object, type = "SnpI")$Color=="Red",,drop=FALSE]
    M.I.Red <- getRed(object)[snpProbesI.Red$AddressB,,drop=FALSE]
    U.I.Red <- getRed(object)[snpProbesI.Red$AddressA,,drop=FALSE]
    M.I.Green <- getGreen(object)[snpProbesI.Green$AddressB,,drop=FALSE]
    U.I.Green <- getGreen(object)[snpProbesI.Green$AddressA,,drop=FALSE]
    
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
              callGeneric(object, ...)
          })

setMethod("mapToGenome", signature(object = "RGChannelSet"),
          function(object, ...) {
              object <- preprocessRaw(object)
              callGeneric(object, ...)
          })

subsetByLoci <- function(rgSet, includeLoci = NULL, excludeLoci = NULL, keepControls = TRUE, keepSnps = TRUE){
    .isRGOrStop(rgSet)
    TypeII <- getProbeInfo(rgSet, type = "II")
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    if(!is.null(includeLoci)) {
        TypeII <- TypeII[TypeII$Name %in% includeLoci, ]
        TypeI.Red <- TypeI.Red[TypeI.Red$Name %in% includeLoci, ]
        TypeI.Green <- TypeI.Green[TypeI.Green$Name %in% includeLoci, ]
    }
    if(!is.null(excludeLoci)) {
        TypeII <- TypeII[! TypeII$Name %in% excludeLoci, ]
        TypeI.Red <- TypeI.Red[! TypeI.Red$Name %in% excludeLoci, ]
        TypeI.Green <- TypeI.Green[! TypeI.Green$Name %in% excludeLoci, ]
    }
    addresses <- c(TypeII$AddressA, 
                   TypeI.Red$AddressA, TypeI.Red$AddressB,
                   TypeI.Green$AddressA, TypeI.Green$AddressB)
    if (keepControls){
        addresses <- c(addresses, getProbeInfo(rgSet, type = "Control")$Address)
    }
    if (keepSnps){
        addresses <- c(addresses,
                       getProbeInfo(rgSet, type = "SnpI")$AddressA, getProbeInfo(rgSet, type = "SnpI")$AddressB,
                       getProbeInfo(rgSet, type = "SnpII")$AddressA)
    }
    indices <- which(rownames(rgSet) %in% addresses)
    rgSet <- rgSet[indices,]
    rgSet
}

setMethod("combine", signature(x = "RGChannelSet", y = "RGChannelSet"),
          function(x, y, ...) {
    colDataFix <- .harmonizeDataFrames(.pDataFix(colData(x)),
                                       .pDataFix(colData(y)))
    colData(x) <- colDataFix$x
    colData(y) <- colDataFix$y
    cbind(x,y)
})
                   
setMethod("coerce", signature(from = "RGChannelSetExtended", to = "RGChannelSet"),
          function(from, to){
    if(nrow(from) > 0 || ncol(from) > 0)
        assays(from) <- assays(from)[c("Green", "Red")]
    class(from) <- "RGChannelSet"
    from
})
