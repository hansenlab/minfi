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
    .show.annotation(object@annotation)
    cat("Number of type I probes:", nrow(object@data[["TypeI"]]), "\n")
    cat("Number of type II probes:", nrow(object@data[["TypeII"]]), "\n")
    cat("Number of control probes:", nrow(object@data[["TypeControl"]]), "\n")
})

IlluminaMethylationManifest <- function(TypeI = new("data.frame"), TypeII = new("data.frame"),
                        TypeControl = new("data.frame"), annotation = "") {
    data <- new.env(parent = emptyenv())
    data[["TypeI"]] <- TypeI
    data[["TypeII"]] <- TypeII
    data[["TypeControl"]] <- TypeControl
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

setClass("IlluminaMethylationAnnotation",
         representation(data = "environment",
                        annotation = "character"))

setValidity("IlluminaMethylationAnnotation", function(object) {
    msg <- NULL
    if (is.null(msg)) TRUE else msg
})

setMethod("show", "IlluminaMethylationAnnotation", function(object) {
    cat("IlluminaMethylationAnnotation object\n")
    .show.annotation(object@annotation)
})

IlluminaMethylationAnnotation <- function(listOfObjects,
                                          annotation = "") {
    data <- new.env(parent = emptyenv())
    stopifnot(all(c("Locations.hg18", "Locations.hg19") %in% names(listOfObjects)))
    for(nam in names(listOfObjects)) {
        cat(nam, "\n")
        assign(nam, listOfObjects[[nam]], envir = data)
    }
    lockEnvironment(data, bindings = TRUE)
    anno <- new("IlluminaMethylationAnnotation",
                annotation = annotation, data = data)
    anno
}

setMethod("getManifest", signature(object = "IlluminaMethylationAnnotation"),
          function(object) {
              maniString <- .getManifestString(object@annotation)
              if(!require(maniString, character.only = TRUE))
                  stop(sprintf("cannot load manifest package %s", maniString))
              get(maniString)
          })

setMethod("getLocations", signature(object = "IlluminaMethylationAnnotation"),
          function(object, genomeBuild = "hg19", mergeManifest = FALSE) {
              locName <- sprintf("Locations.%s", genomeBuild)
              if(! locName %in% ls(object@data))
                  stop(sprintf("genomeBuild '%s' not in object", genomeBuild))
              locations <- get(locName, envir = object@data)
              locations[locations[, "chr"] == "chr", "chr"] <- "unmapped"
              wh.unmap <- which(locations[, "chr"] == "unmapped")
              if(length(wh.unmap) > 0)
                  locations[wh.unmap, "pos"] <- seq(from = 1, by = 2, length.out = length(wh.unmap))
              if("strand" %in% colnames(locations))
                  strand <- locations[, "strand"]
              else
                  strand <- rep("*", nrow(locations))
              gr <- GRanges(seqnames = locations[, "chr"],
                            strand = strand,
                            ranges = IRanges(start = locations[, "pos"], width = 1))
              genome(gr) <- genomeBuild
              names(gr) <- rownames(locations)
              if(mergeManifest) {
                  typeI <- getProbeInfo(object, type = "I")
                  typeI$Type <- "I"
                  typeII <- getProbeInfo(object, type = "II")
                  typeII$Type <- "II"
                  names(typeII)[names(typeII) == "ProbeSeq"] <- "ProbeSeqA"
                  names(typeII)[names(typeII) == "Address"] <- "AddressA"
                  typeII$AddressB <- rep(NA_character_, nrow(typeII))
                  typeII$Color <- rep(NA_character_, nrow(typeII))
                  typeII$NextBase <- rep(NA_character_, nrow(typeII))
                  typeII$ProbeSeqB <- rep(NA_character_, nrow(typeII))
                  manifest <- as(rbind(typeI, typeII[, names(typeI)]), "DataFrame")
                  rownames(manifest) <- manifest$Name
                  manifest[, "ProbeSeqA"] <- DNAStringSet(manifest[, "ProbeSeqA"])
                  manifest[, "ProbeSeqB"] <- DNAStringSet(manifest[, "ProbeSeqB"])
                  manifest <- manifest[names(gr),]
                  elementMetadata(gr) <- manifest
              }
              gr
          })
