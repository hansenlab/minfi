getManifest <- function(object) {
    name <- paste(annotation(object), "manifest", sep = "")
    get(name)
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

getRed <- function(object) {
    if(!is(object, "RGChannelSet"))
        stop("'object' needs to be a 'RGChannelSet'")
    assayDataElement(object, "Red")
}

getGreen <- function(object) {
    if(!is(object, "RGChannelSet"))
        stop("'object' needs to be a 'RGChannelSet'")
    assayDataElement(object, "Green")
}

getMeth <- function(object) {
    if(!is(object, "MethylSet"))
        stop("'object' needs to be an 'MethylSet'")
    assayDataElement(object, "Meth")
}

getUnmeth <- function(object) {
    if(!is(object, "MethylSet"))
        stop("'object' needs to be an 'MethylSet'")
    assayDataElement(object, "Unmeth")
}

getBeta <- function(object, type = c("minfi", "Illumina")) {
    if(!is(object, "MethylSet") & !is(object, "RGChannelSet"))
        stop("'object' needs to be an 'MethylSet' or a 'RGChannelSet'")
    type <- match.arg(type)
    if(is(object, "RGChannelSet"))
        object <- preprocessRaw(object)
    Meth <- pmax(assayDataElement(object, "Meth"), 0)
    Unmeth <- pmax(assayDataElement(object, "Unmeth"), 0)
    if(type == "minfi")
        return(Meth / (Meth + Unmeth))
    if(type == "Illumina")
        return(Meth / (Meth + Unmeth + 100))
}

getM <- function (object, type = c("minfi", "Illumina")) {
    if(!is(object, "MethylSet") & !is(object, "RGChannelSet"))
        stop("'object' needs to be an 'MethylSet' or a 'RGChannelSet'")
    type <- match.arg(type)
    beta <- getBeta(object, type = type)
    beta <- pmin(beta, 0.999)
    beta <- pmax(beta, 0.001)
    logit2(beta)
}

