getProbeInfo <- function(object, type = c("I", "II", "Control", "I-Green", "I-Red")) {
    type <- match.arg(type)
    if(type %in% c("I", "II", "Control"))
        return(getManifest(object)@data[[paste("Type", type, sep = "")]])
    typeI <- getManifest(object)@data[["TypeI"]]
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

