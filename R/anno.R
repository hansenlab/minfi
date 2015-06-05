setClass("IlluminaMethylationAnnotation",
         representation(data = "environment",
                        annotation = "character",
                        defaults = "character"
                        ))

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
    .show.availableAnnotation(object)
})

IlluminaMethylationAnnotation <- function(listOfObjects, annotation = "", defaults = "") {
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
    stopifnot(all(listOfObjects[["Locations"]]$chr %in% .seqnames.order.all))
    available <- .availableAnnotation(listOfObjects)
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

getAnnotationObject <- function(object) {
    if(is(object, "MethylSet") || is(object, "RatioSet") ||
       is(object, "GenomicMethylSet") || is(object, "GenomicRatioSet") ||
       is(object, "RGChannelSet"))
        object <- .getAnnotationString(object@annotation)
    if(is.character(object)) {
        if(!require(object, character.only = TRUE))
            stop(sprintf("cannot load annotation package %s", object))
        object <- get(object)
    }
    if(!is(object, "IlluminaMethylationAnnotation"))
        stop("Could not locate annotation object for 'object' of class", class(object))
    object
}

.availableAnnotation <- function(object) {
    object <- getAnnotationObject(object)
    allAnnoNames <- ls(object@data)
    annoClasses <- sub("\\..*", "", allAnnoNames)
    annoClassesChoices <- sub(".*\\.", "", allAnnoNames)
    annoClassesChoices[grep("\\.", allAnnoNames, invert = TRUE)] <- ""
    annoClassesChoices <- split(annoClassesChoices, annoClasses)
    annoClasses <- unique(annoClasses)
    defaults <- object@defaults
    out <- list(names = allAnnoNames, classes = annoClasses,
                classChoices = annoClassesChoices, defaults = defaults)
    out
}

getAnnotation <- function(object, what = "everything", lociNames = NULL,
                          orderByLocation = FALSE, dropNonMapping = FALSE) {
    ## processing of arguments and check
    annoObject <- getAnnotationObject(object)
    available <- .availableAnnotation(annoObject)
    if("everything" %in% what)
        what <- available$defaults
    if(!(all(what %in% available$names)))
        stop("the value of argument 'what' is not part of the annotation package or 'everything'")
    if(any(sapply(available$classes, function(cl) {length(grep(cl, what))} ) > 1))
        stop("only on choice per allowable class")
    if(!any(grepl("^Locations", what)) && (orderByLocation || dropNonMapping))
        stop("To use 'orderbyLocation' or 'dropNonMapping' specify Locations as part of 'what'")
    ## FIXME: ensure same order always
    ## Old code for inspiration
    ## bestOrder <- c("Locations", "Manifest", "IlluminaSNPs", "Annotation")
    ## what <- bestOrder[bestOrder %in% what]
    
    out <- do.call(cbind, lapply(what, function(wh) {
        get(wh, envir = annoObject@data)
    }))

    if(!is.null(lociNames)) {
        lociNames <- lociNames[lociNames %in% rownames(out)]
    }
    if(is(object, "MethylSet") || is(object, "RatioSet") ||
       is(object, "GenomicMethylSet") || is(object, "GenomicRatioSet")) {
        fNames <- featureNames(object)
        if(is.null(lociNames))
            lociNames <- fNames[fNames %in% rownames(out)]
        else
            lociNames <- fNames[fNames %in% lociNames]
    }
    if(!is.null(lociNames))
        out <- out[match(lociNames, rownames(out)),]
    if(dropNonMapping) {
        seqOrder <- .seqnames.order
        wh <- which(out$chr %in% seqOrder)
        out <- out[wh,]
    } else {
        seqOrder <- .seqnames.order.all
    }
    if(orderByLocation) {
        ## we have a data.frame with chr, pos.  We want to sort accoridng to chr nested with pos, but with a custom
        ## ordering on chr.
        ## we will make chr into an integer
        chrInteger <- as.integer(factor(out$chr, levels = seqOrder))
        od <- order(chrInteger, out$pos)
        out <- out[od,]
        ## ## FIXME: very memory intensive
        ## sp <- split(out, out$chr)
        ## sp <- sp[seqOrder[seqOrder %in% names(sp)]]
        ## out <- do.call(rbind, lapply(sp, function(df) {
        ##     od <- order(df$pos)
        ##     df[od,]
        ## }))
    }
    out
}

getLocations <- function(object, mergeManifest = FALSE,
                         orderByLocation = FALSE, lociNames = NULL) {
    if(mergeManifest)
        what <- c("Locations", "Manifest")
    else
        what <- "Locations"
    locs <- getAnnotation(object, what = what, dropNonMapping = TRUE,
                          orderByLocation = orderByLocation, lociNames = lociNames)
    gr <- GRanges(seqnames = locs$chr,
                  ranges = IRanges(start = locs$pos, width = 1))
    seqlevels(gr) <- .seqnames.order[.seqnames.order %in% seqlevels(gr)]
    names(gr) <- rownames(locs)
    if(mergeManifest) {
        if("strand" %in% names(locs))
            names(locs)[names(locs) == "strand"] <- "assayStrand"
        mcols(gr) <- locs[, ! names(locs) %in% c("chr", "pos", "strand")]
    }
    genome(gr) <- unname(getAnnotationObject(object)@annotation["genomeBuild"])
    gr
}

.getIslandAnnotation <- function(object, islandAnno = NULL) {
    av <- .availableAnnotation(object)
    if(is.null(islandAnno)) {
        islandAnno <- grep("^Islands\\.", getAnnotationObject(object)@defaults, value = TRUE)
    } else {
        islandAnno <- sub("^Islands\\.", "", islandAnno)
        if(! islandAnno %in% av$annoClassesChoices) {
            stop(sprintf("islandAnno '%s' is not part of the annotation", islandAnno))
        } else {
            islandAnno <- sprintf("Islands.%s", islandAnno)
        }
    }
    getAnnotation(object, what = islandAnno)
}
    
getIslandStatus <- function(object, islandAnno = NULL) {
    regionType <- .getIslandAnnotation(object, islandAnno = islandAnno)$Relation_to_Island
    regionType <- sub("^[SN]_", "", regionType)
    regionType
}

getProbeType <- function(object, withColor = FALSE) {
    if(withColor) {
        probeType <- paste0(getAnnotation(object, what = "Manifest")$Type,
                           getAnnotation(object, what = "Manifest")$Color)
    } else {
        probeType <- getAnnotation(object, what = "Manifest")$Type
    }
    probeType
}

getSnpInfo <- function(object, snpAnno = NULL) {
    av <- .availableAnnotation(object)
    if(is.null(snpAnno)) {
        snpAnno <- grep("^SNPs\\.", getAnnotationObject(object)@defaults, value = TRUE)
    } else {
        snpAnno <- sub("^SNPs\\.", "", snpAnno)
        if(! snpAnno %in% av$annoClassesChoices) {
            stop(sprintf("snpAnno '%s' is not part of the annotation", snpAnno))
        } else {
            snpAnno <- sprintf("SNPs.%s", snpAnno)
        }
    }
    snps <- getAnnotation(object, what = snpAnno)
    snps
}

addSnpInfo <- function(object, snpAnno = NULL) {
    .isGenomic(object)
    snps <- getSnpInfo(object = object, snpAnno = snpAnno)
    mcolsNames <- names(mcols(granges(object)))
    if(any(mcolsNames %in% names(snps)))
        cat("Replacing existing snp information\n")
    mcols(object@rowRanges) <- cbind(mcols(granges(object)), snps)
    object
}

.doSnpOverlap <- function(map, grSnp) {
    stopifnot(is(map, "GRanges"))
    stopifnot(is(grSnp, "GRanges"))
    stopifnot(all(c("SBE", "probeStart", "probeEnd") %in% names(mcols(map))))
    cat("removing Snps with width != 1\n")
    grSnp <- grSnp[width(grSnp) == 1]
    cpgGR <- GRanges(seqnames(map), IRanges(start(map), width=2))
    ooCpG <- findOverlaps(cpgGR, grSnp)
    sbeGR <- GRanges(seqnames(map), IRanges(map$SBE, map$SBE))
    ooSbe <- findOverlaps(sbeGR, grSnp)
    ## just match to whole probe then drop CpG overlaps
    probeGR <- GRanges(seqnames(map), IRanges(map$probeStart, map$probeEnd))
    ooProbe <- findOverlaps(probeGR, grSnp, ignore.strand=TRUE)
    ooProbe <- ooProbe[-which(queryHits(ooProbe) %in% queryHits(ooCpG))]
    snpAnno <- DataFrame(matrix(nrow = length(map), ncol = 6))
    colnames(snpAnno) = c("Probe_rs" , "Probe_maf", "CpG_rs",
            "CpG_maf" ,  "SBE_rs" ,   "SBE_maf")
    rownames(snpAnno) <- names(map)
    snpAnno$Probe_rs[queryHits(ooProbe)] <- names(grSnp)[subjectHits(ooProbe)]
    snpAnno$Probe_maf[queryHits(ooProbe)] <- grSnp$MAF[subjectHits(ooProbe)]
    snpAnno$CpG_rs[queryHits(ooCpG)] <- names(grSnp)[subjectHits(ooCpG)]
    snpAnno$CpG_maf[queryHits(ooCpG)] <- grSnp$MAF[subjectHits(ooCpG)]
    snpAnno$SBE_rs[queryHits(ooSbe)] <- names(grSnp)[subjectHits(ooSbe)]
    snpAnno$SBE_maf[queryHits(ooSbe)] <- grSnp$MAF[subjectHits(ooSbe)]
    snpAnno
}

dropLociWithSnps <- function(object, snps = c("CpG", "SBE"), maf = 0, snpAnno = NULL){
    .isGenomic(object)
    maf_cols <- paste0(snps, "_maf")
    snpDF  <- getSnpInfo(object, snpAnno = snpAnno)
    choices <- c("Probe_maf", "CpG_maf", "SBE_maf")
    if(!all(choices %in% colnames(snpDF)))
        stop("The specificed 'snpAnno' is not supported by this function")
    if (sum(!(maf_cols %in% choices))>0){
        stop("snps vector argument must be a combination of  \"Probe\", \"CpG\" and \"SBE\"")
    }
    if (!is.numeric(maf) || maf<0 || maf>1){
        stop("maf argument must be a numeric value between 0 and 1")
    }
    wh <- Reduce(union, lapply(maf_cols, function(xx) {
        which(snpDF[, xx] >= maf)
    }))
    wh <- sort(wh)
    if(length(wh) == 0)
        return(object)
    object[-wh,]
}
