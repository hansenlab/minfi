getManifest <- function(object) {
    maniString <- .getManifestString(object@annotation)
    if(!require(maniString, character.only = TRUE))
        stop(sprintf("cannot load manifest package %s", maniString))
    get(maniString)
}

getProbeData <- function(object) {
    if(is(object, "IlluminaMethylationManifest"))
        return(object@data)
    if(is(object, "RGChannelSet"))
        return(getManifest(object)@data)
    if(is(object, "MethylSet"))
        return(getManifest(object)@data)
    stop("cannot handle 'object'")
}

getProbeInfo <- function(object, type = c("I", "II", "Control", "I-Green", "I-Red")) {
    type <- match.arg(type)
    if(type %in% c("I", "II", "Control"))
        return(getProbeData(object)[[paste("Type", type, sep = "")]])
    typeI <- getProbeData(object)[["TypeI"]]
    if(type == "I-Green")
        return(typeI[typeI$Color == "Grn",])
    if(type == "I-Red")
        return(typeI[typeI$Color == "Red",])
}

getManifestInfo <- function(object, type = c("nLoci", "locusNames")) {
    type <- match.arg(type)
    switch(type,
           "nLoci" = {
               nrow(getProbeInfo(object, type = "I")) +
                   nrow(getProbeInfo(object, type = "II"))
           },
           "locusNames" = {
               c(getProbeInfo(object, type = "I")$Name,
                 getProbeInfo(object, type = "II")$Name)
           })
}

getControlAddress <- function(object, controlType = c("NORM_A", "NORM_C", "NORM_G", "NORM_T")) {
    ctrls <- getProbeInfo(object, type = "Control")
    ctrls[ctrls$Type %in% controlType, "Address"]
}

