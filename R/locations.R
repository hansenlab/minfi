setMethod("mapToGenome", signature(object = "RGChannelSet"),
          function(object, ...) {
              object <- preprocessRaw(object)
              callGeneric(object, ...)
          })

setMethod("mapToGenome", signature(object = "MethylSet"),
          function(object, genomeBuild = c("hg19", "hg18"),
                   drop = TRUE, mergeManifest = FALSE) {
              genomeBuild <- match.arg(genomeBuild)
              object <- dropMethylationLoci(object, dropRS = TRUE, dropCH = FALSE)
              gr <- getLocations(object, genomeBuild = genomeBuild, drop = drop,
                                 mergeManifest = mergeManifest)
              gr <- sort(gr)
              object <- object[names(gr),]
              GenomicMethylSet(gr = gr, Meth = getMeth(object),
                               Unmeth = getUnmeth(object),
                               pData = pData(object),
                               preprocessMethod = preprocessMethod(object),
                               annotation = annotation(object))
          })

getAnnotation <- function(object, genomeBuild = c("hg19", "hg18"), what = "everything",
                          returnAs = c("data.frame", "GRanges"),
                          orderByLocation = FALSE, drop = FALSE) {
    ## processing of arguments and check
    returnAs <- match.arg(returnAs)
    genomeBuild <- match.arg(genomeBuild)
    annoString <- .getAnnotationString(object@annotation)
    if(!require(annoString, character.only = TRUE))
        stop(sprintf("cannot load annotation package %s", annoString))
    annotation <- get(annoString)
    possible <- unique(sub("\\.hg.*$", "", ls(annotation@data)))
    if("everything" %in% what)
        what <- possible
    if(!(all(what %in% possible)))
        stop("the value of argument 'what' is not part of the annotation package or 'everything'")
    if(returnAs == "GRanges" && ! ("Locations" %in% what))
        what <- c("Locations", what)
    locName <- paste("Locations", genomeBuild, sep = ".")
    bestOrder <- c("Locations", "Manifest", "IlluminaSNPs", "Annotation")
    what <- bestOrder[bestOrder %in% what]
    if("Locations" %in% what)
        what[what == "Locations"] <- locName
    if(!is(object, "MethylSet"))
        stop(sprintf("'objects' is of class '%s' which is not yet supported", class(object)))
    
    ## We need locs to filter out loci using the function arguments 'order' and 'drop'
    locs <- get(locName, annotation@data)
    locs <- locs[featureNames(object),]
    locs$Name <- rownames(locs)
    if(orderByLocation == TRUE) {
        if(drop) {
            chrs <- paste("chr", c(1:22, "X", "Y"), sep = "")
        } else {
            chrs <- paste("chr", c(1:22, "X", "Y", ""), sep = "")
        }
        sp <- split(locs, locs$chr)[chrs]
        locs <- do.call(rbind, lapply(sp, function(df) {
            od <- order(df$pos)
            df[od,]
        }))
        rownames(locs) <- locs$Name
        locs$Name <- NULL
    } else {
        if(drop)
            locs <- locs[locs$chr != "chr",]
    }
    out <- do.call(cbind, lapply(what, function(wh) {
        get(wh, envir = annotation@data)
    }))
    out <- out[rownames(locs),]
    if(returnAs == "data.frame") {
        return(out)
    }
    if(returnAs == "GRanges")
        return(GRanges(seqnames = out$chr, ranges = IRanges(start = out$pos, width = 2),
                       out[, !names(out) %in% c("chr", "pos")]))
}

orderByLocation <- function(x, what = "everything", genomeBuild = "hg19",
                              returnAs = "data.frame", drop = TRUE) {
    anno <- getAnnotation(object = x, what = what,
                          genomeBuild = genomeBuild, returnAs = returnAs,
                          orderByLocation = TRUE, drop = drop)
    out <- list(meth = getMeth(x)[rownames(anno),],
                unmeth = getUnmeth(x)[rownames(anno),],
                locs = anno[, c("chr", "pos")],
                everything = anno)
    return(out)
}

