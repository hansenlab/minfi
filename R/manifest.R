setClass("IlluminaMethylationManifest",
         representation(data = "environment",
                        annotation = "character"))

setValidity("IlluminaMethylationManifest", function(object) {
    msg <- NULL
    if(! "TypeI" %in% ls(object@data) || !is(object@data[["TypeI"]], "DataFrame"))
        msg <- paste(msg, "'data' slot must contain a DataFrame with TypeI probes", sep = "\n")
    if(! "TypeII" %in% ls(object@data) || !is(object@data[["TypeII"]], "DataFrame"))
        msg <- paste(msg, "'data' slot must contain a DataFrame with TypeII probes", sep = "\n")
    if(! "TypeControl" %in% ls(object@data) || !is(object@data[["TypeControl"]], "DataFrame"))
        msg <- paste(msg, "'data' slot must contain a DataFrame with Control probes", sep = "\n")
    if(! "TypeSnpI" %in% ls(object@data) || !is(object@data[["TypeSnpI"]], "DataFrame"))
        msg <- paste(msg, "'data' slot must contain a DataFrame with Snp I probes", sep = "\n")
    if(! "TypeSnpII" %in% ls(object@data) || !is(object@data[["TypeSnpII"]], "DataFrame"))
        msg <- paste(msg, "'data' slot must contain a DataFrame with Snp II probes", sep = "\n")

    ## Check Names
    if(! all(c("Name", "AddressA", "AddressB", "Color", "nCpG") %in% colnames(object@data[["TypeI"]])))
        msg <- paste(msg, "'TypeI' has to have column names 'Name', 'AddressA', 'AddressB', 'Color', 'nCpG'")
    if(!is.character(object@data[["TypeI"]]$Name) ||
       !is.character(object@data[["TypeI"]]$AddressA) ||
       !is.character(object@data[["TypeI"]]$AddressB) ||
       !is.character(object@data[["TypeI"]]$Color) ||
       !is.integer(object@data[["TypeI"]]$nCpG))
        msg <- paste(msg, "'TypeI' columns has wrong classes")
    
    if(! all(c("Name", "AddressA", "nCpG") %in% colnames(object@data[["TypeII"]])))
        msg <- paste(msg, "'TypeII' has to have column names 'Name', 'AddressA', 'nCpG'")
    if(!is.character(object@data[["TypeII"]]$Name) ||
       !is.character(object@data[["TypeII"]]$AddressA) ||
       !is.integer(object@data[["TypeII"]]$nCpG))
        msg <- paste(msg, "'TypeII' columns has wrong classes")
    
    if(! all(c("Type", "Address") %in% colnames(object@data[["TypeControl"]])))
        msg <- paste(msg, "'TypeControl' has to have column names 'Type', 'Address'")
    if(!is.character(object@data[["TypeControl"]]$Type) ||
       !is.character(object@data[["TypeControl"]]$Address))
        msg <- paste(msg, "'TypeControl' columns has wrong classes")
    
    if(! all(c("Name", "AddressA", "AddressB", "Color") %in% colnames(object@data[["TypeSnpI"]])))
        msg <- paste(msg, "'TypeSnpI' has to have column names 'Name', 'AddressA', 'AddressB', 'Color'")
    if(!is.character(object@data[["TypeSnpI"]]$Name) ||
       !is.character(object@data[["TypeSnpI"]]$AddressA) ||
       !is.character(object@data[["TypeSnpI"]]$AddressB) ||
       !is.character(object@data[["TypeSnpI"]]$Color))
        msg <- paste(msg, "'TypeSnpI' columns has wrong classes")

    if(! all(c("Name", "AddressA") %in% colnames(object@data[["TypeSnpII"]])))
        msg <- paste(msg, "'TypeSnpII' has to have column names 'Name', 'AddressA'")
    if(!is.character(object@data[["TypeSnpII"]]$Name) ||
       !is.character(object@data[["TypeSnpII"]]$AddressA))
        msg <- paste(msg, "'TypeSnpII' columns has wrong classes")
    
    if (is.null(msg)) TRUE else msg
})

setMethod("show", "IlluminaMethylationManifest", function(object) {
    cat("IlluminaMethylationManifest object\n")
    .show.annotation(object@annotation)
    cat("Number of type I probes:", nrow(object@data[["TypeI"]]), "\n")
    cat("Number of type II probes:", nrow(object@data[["TypeII"]]), "\n")
    cat("Number of control probes:", nrow(object@data[["TypeControl"]]), "\n")
    cat("Number of SNP type I probes:", nrow(object@data[["TypeSnpI"]]), "\n")
    cat("Number of SNP type II probes:", nrow(object@data[["TypeSnpII"]]), "\n")
})

IlluminaMethylationManifest <- function(TypeI = new("DataFrame"),
                                        TypeII = new("DataFrame"),
                                        TypeControl = new("DataFrame"),
                                        TypeSnpI = new("DataFrame"),
                                        TypeSnpII = new("DataFrame"),
                                        annotation = "") {
    data <- new.env(parent = emptyenv())
    data[["TypeI"]] <- TypeI
    data[["TypeII"]] <- TypeII
    data[["TypeControl"]] <- TypeControl
    data[["TypeSnpI"]] <- TypeSnpI
    data[["TypeSnpII"]] <- TypeSnpII
    lockEnvironment(data, bindings = TRUE)
    manifest <- new("IlluminaMethylationManifest", annotation = annotation, data = data)
    manifest
}

setMethod("getManifest", signature(object = "IlluminaMethylationManifest"),
          function(object) {
              object
          })

setMethod("getManifest", signature(object = "character"),
          function(object) {
              maniString <- .getManifestString(object)
              if(!require(maniString, character.only = TRUE))
                  stop(sprintf("cannot load manifest package %s", maniString))
              get(maniString)
          })

getProbeInfo <- function(object, type = c("I", "II", "Control", "I-Green", "I-Red", "SnpI", "SnpII")) {
    type <- match.arg(type)
    if(type %in% c("I", "II", "Control", "SnpI", "SnpII"))
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

getControlAddress <- function(object, controlType = c("NORM_A", "NORM_C", "NORM_G", "NORM_T"),
                              asList = FALSE) {
    if(asList) {
        ctrls <- getProbeInfo(object, type = "Control")
        out <- split(ctrls$Address, ctrls$Type)
        out <- out[names(out) %in% controlType]
    } else {
        ctrls <- getProbeInfo(object, type = "Control")
        out <- ctrls[ctrls$Type %in% controlType, "Address"]
    }
    out
}

getControlTypes <- function(object) {
    ctrls <- getProbeInfo(object, type = "Control")
    table(ctrls[, ""])
}

getProbePositionsDetailed <- function(map) {
    ## map is GR with metadata columns strand and type
    stopifnot(is(map, "GRanges"))
    stopifnot(c("Strand", "Type") %in% names(elementMetadata(map)))

    probeStart <- rep(NA, length(map))
    wh.II.F <- which(map$Type=="II" & map$Strand=="+")
    wh.II.R <- which(map$Type=="II" & map$Strand=="-")
    wh.I.F <- which(map$Type=="I" & map$Strand=="+")
    wh.I.R <- which(map$Type=="I" & map$Strand=="-")
    
    probeStart[wh.II.F] <- start(map)[wh.II.F]
    probeStart[wh.II.R] <- start(map)[wh.II.R] - 50
    probeStart[wh.I.F] <- start(map)[wh.I.F] - 1
    probeStart[wh.I.R] <- start(map)[wh.I.R] - 49
    map$probeStart <- probeStart

    probeEnd <- rep(NA, length(map))
    probeEnd[wh.II.F] <- start(map)[wh.II.F] + 50 
    probeEnd[wh.II.R] <- start(map)[wh.II.R] 
    probeEnd[wh.I.F] <- start(map)[wh.I.F] + 49
    probeEnd[wh.I.R] <- start(map)[wh.I.R] + 1
    map$probeEnd <- probeEnd
    
    sbe <- rep(NA, length(map))
    sbe[wh.II.F] <- start(map)[wh.II.F] 
    sbe[wh.II.R] <- start(map)[wh.II.R] + 1
    sbe[wh.I.F] <- start(map)[wh.I.F] - 1
    sbe[wh.I.R] <- start(map)[wh.I.R] + 2
    map$SBE <- sbe

    map
}
    
