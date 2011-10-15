setClass("RGChannelSet", contains = "eSet")

setValidity("RGChannelSet", function(object) {
    msg <- validMsg(NULL, assayDataValidMembers(assayData(object),
                                                c("Red", "Green")))
    if (is.null(msg)) TRUE else msg
})

setMethod("initialize", "RGChannelSet",
          function(.Object, Red = new("matrix"), Green = new("matrix"), ...) {
              callNextMethod(.Object, Red = Red, Green = Green, ...)
})

setClass("RGChannelSetExtended", contains = "RGChannelSet")

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

setClass("MethylSet",
         representation(preprocessMethod = "character"),
         contains = "eSet")

setValidity("MethylSet", function(object) {
    msg <- validMsg(NULL, assayDataValidMembers(assayData(object),
                                                c("Meth", "Unmeth")))
    if (is.null(msg)) TRUE else msg
})

setMethod("initialize", "MethylSet",
          function(.Object, Meth = new("matrix"), Unmeth = new("matrix"), ...){
              callNextMethod(.Object, Meth = Meth, Unmeth = Unmeth, ...)
})

setMethod("show", "MethylSet", function(object) {
    callNextMethod(object)
    preprocess <- object@preprocessMethod
    if(length(preprocess) == 0)
        preprocess <- c("unknown", "unknown", "unknown")
    cat("Preprocessing", preprocess[1],
        "minfi version", preprocess[2],
        "Manifest version", preprocess[3], fill = TRUE)
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

getManifest <- function(object) {
    name <- paste(annotation(object), "manifest", sep = "")
    get(name)
}

getProbeData <- function(object) {
    if(is(object, "IlluminaMethylationManifest"))
        return(object@data)
    if(is(object, "RGChannelSet"))
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

