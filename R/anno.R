setClass("IlluminaMethylationAnnotation",
         representation(data = "environment",
                        annotation = "character",
                        defaults = "character"))

setValidity("IlluminaMethylationAnnotation", function(object) {
    msg <- NULL
    if(!(all(sapply(object@data, class) == "DataFrame")))
        msg <- paste(msg,
                     "All objects in 'objects@data' has to be of class 'DataFrame'", sep = "\n")
    if (is.null(msg)) TRUE else msg
})

setMethod("show", "IlluminaMethylationAnnotation", function(object) {
    cat("IlluminaMethylationAnnotation object\n")
    .show.annotation(object@annotation)
})

IlluminaMethylationAnnotation <- function(listOfObjects, annotation = "",
                                          defaults = "") {
    stopifnot(annotation != "")
    stopifnot(all(c("array", "annotation", "genomeBuild") %in% names(annotation)))
    stopifnot(all(c("Manifest", "Locations") %in% names(listOfObjects)))
    Manifest <- listOfObjects[["Manifest"]]
    stopifnot(setequal(names(Manifest),
                       c("Name", "AddressA", "AddressB", "ProbeSeqA",
                         "ProbeSeqB", "Type", "NextBase", "Color")))
    stopifnot(all(sapply(listOfObjects, class) %in% c("DataFrame", "data.frame")))
    stopifnot(all(nrow(Manifest) == sapply(listOfObjects, nrow)))
    stopifnot(all(sapply(listOfObjects, function(obj) {
        all(rownames(obj) == rownames(Manifest))
    })))
    stopifnot(all(c("chr", "pos") %in% names(listOfObjects[["Locations"]])))
    stopifnot(all(Locations$chr %in% minfi:::.seqnames.order.all))
    available <- minfi:::availableAnnotation(listOfObjects, constructorCheck = TRUE)
    stopifnot(all(defaults %in% names(listOfObjects)))
    stopifnot(!anyDuplicated(sub("\\..*", "", defaults)))
    ## FIXME: Check column names of any Islands object
    
    ## Instantiating
    data <- new.env(parent = emptyenv())
    for(nam in names(listOfObjects)) {
        cat(nam, "\n")
        assign(nam, as(listOfObjects[[nam]], "DataFrame"), envir = data)
    }
    lockEnvironment(data, bindings = TRUE)
    anno <- new("IlluminaMethylationAnnotation",
                annotation = annotation, data = data, defaults = defaults)
    anno
}

setMethod("getManifest", signature(object = "IlluminaMethylationAnnotation"),
          function(object) {
              maniString <- .getManifestString(object@annotation)
              if(!require(maniString, character.only = TRUE))
                  stop(sprintf("cannot load manifest package %s", maniString))
              get(maniString)
          })

setMethod("getLocations", signature(object = "character"),
          function(object, mergeManifest = FALSE) {
              callNextMethod(get(object), mergeManifest = mergeManifest)
          })

setMethod("getLocations", signature(object = "IlluminaMethylationAnnotation"),
          function(object, mergeManifest = FALSE) {
              locName <- "Locations"
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
              ## FIXME
              ## genome(gr) <- genomeBuild
              names(gr) <- rownames(locations)
              if(mergeManifest) {
                  typeI <- rbind(getProbeInfo(object, type = "I"),
                                 getProbeInfo(object, type = "SnpI"))
                  typeI$Type <- "I"
                  typeII <- rbind(getProbeInfo(object, type = "II"),
                                  getProbeInfo(object, type = "SnpII"))
                  typeII$Type <- "II"
                  typeII$AddressB <- rep(NA_character_, nrow(typeII))
                  typeII$Color <- rep(NA_character_, nrow(typeII))
                  typeII$NextBase <- DNAStringSet(character(nrow(typeII)))
                  typeII$ProbeSeqB <- DNAStringSet(character(nrow(typeII)))
                  stopifnot(setequal(names(typeI), names(typeII)))
                  manifest <- rbind(typeI, typeII[, names(typeI)])
                  rownames(manifest) <- manifest$Name
                  manifest <- manifest[names(gr),]
                  elementMetadata(gr) <- manifest
              }
              gr
          })

getIslandStatus <- function(object, islandType = "UCSC") {
    regionType <- minfi:::getAnnotation(object, what = sprintf("Islands.%s", islandType))$Relation_to_Island
    regionType <- sub("^[SN]_", "", regionType)
    regionType
}

doSnpOverlap <- function(map, grSnp) {
    stopifnot(is(map, "GRanges"))
    stopifnot(is(grSnp, "GRanges"))
    stopifnot(all(c("SBE", "probeStart", "probeEnd") %in% names(elementMetadata(map))))
    cpgGR <- GRanges(seqnames(map), IRanges(start(map), width=2))
    ooCpG <- findOverlaps(cpgGR, snpGR)
    sbeGR <- GRanges(seqnames(map), IRanges(map$SBE, map$SBE))
    ooSbe <- findOverlaps(sbeGR, snpGR)
    ## just match to whole probe then drop CpG overlaps
    probeGR <- GRanges(seqnames(map), IRanges(map$probeStart, map$probeEnd))
    ooProbe <- findOverlaps(probeGR, snpGR, ignore.strand=TRUE)
    ooProbe <- ooProbe[-which(queryHits(ooProbe) %in% queryHits(ooCpG))]
    snpAnno <- DataFrame(matrix(nr = length(map), nc = 6))
    colnames(snpAnno) = c("Probe_rs" , "Probe_maf", "CpG_rs",
            "CpG_maf" ,  "SBE_rs" ,   "SBE_maf")
    rownames(snpAnno) <- names(map)
    snpAnno$Probe_rs[queryHits(ooProbe)] <- names(snpGR)[subjectHits(ooProbe)]
    snpAnno$Probe_maf[queryHits(ooProbe)] <- snpGR$MAF[subjectHits(ooProbe)]
    snpAnno$CpG_rs[queryHits(ooCpG)] <- names(snpGR)[subjectHits(ooCpG)]
    snpAnno$CpG_maf[queryHits(ooCpG)] <- snpGR$MAF[subjectHits(ooCpG)]
    snpAnno$SBE_rs[queryHits(ooSbe)] <- names(snpGR)[subjectHits(ooSbe)]
    snpAnno$SBE_maf[queryHits(ooSbe)] <- snpGR$MAF[subjectHits(ooSbe)]
    snpAnno
}
