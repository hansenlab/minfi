## getAnnotation <- function(object, genomeBuild = c("hg19", "hg18"), what = "everything",
##                           returnAs = c("data.frame", "GRanges"),
##                           orderByLocation = FALSE, drop = FALSE) {
##     ## processing of arguments and check
##     returnAs <- match.arg(returnAs)
##     genomeBuild <- match.arg(genomeBuild)
##     annoString <- .getAnnotationString(object@annotation)
##     if(!require(annoString, character.only = TRUE))
##         stop(sprintf("cannot load annotation package %s", annoString))
##     annotation <- get(annoString)
##     possible <- unique(sub("\\.hg.*$", "", ls(annotation@data)))
##     if("everything" %in% what)
##         what <- possible
##     if(!(all(what %in% possible)))
##         stop("the value of argument 'what' is not part of the annotation package or 'everything'")
##     if(returnAs == "GRanges" && ! ("Locations" %in% what))
##         what <- c("Locations", what)
##     locName <- paste("Locations", genomeBuild, sep = ".")
##     bestOrder <- c("Locations", "Manifest", "IlluminaSNPs", "Annotation")
##     what <- bestOrder[bestOrder %in% what]
##     if("Locations" %in% what)
##         what[what == "Locations"] <- locName
##     ## if(!is(object, "MethylSet"))
##     ##     stop(sprintf("'objects' is of class '%s' which is not yet supported", class(object)))
    
##     ## We need locs to filter out loci using the function arguments 'order' and 'drop'
##     locs <- get(locName, annotation@data)
##     locs <- locs[featureNames(object),]
##     locs$Name <- rownames(locs)
##     if(orderByLocation == TRUE) {
##         if(drop) {
##             chrs <- paste("chr", c(1:22, "X", "Y"), sep = "")
##         } else {
##             chrs <- paste("chr", c(1:22, "X", "Y", ""), sep = "")
##         }
##         sp <- split(locs, locs$chr)[chrs]
##         locs <- do.call(rbind, lapply(sp, function(df) {
##             od <- order(df$pos)
##             df[od,]
##         }))
##         rownames(locs) <- locs$Name
##         locs$Name <- NULL
##     } else {
##         if(drop)
##             locs <- locs[locs$chr != "chr",]
##     }
##     out <- do.call(cbind, lapply(what, function(wh) {
##         get(wh, envir = annotation@data)
##     }))
##     out <- out[rownames(locs),]
##     if(returnAs == "data.frame") {
##         return(out)
##     }
##     if(returnAs == "GRanges")
##         return(GRanges(seqnames = out$chr, ranges = IRanges(start = out$pos, width = 2),
##                        out[, !names(out) %in% c("chr", "pos")]))
## }


availableAnnotation <- function(object) {
    object <- minfi:::.getAnnotationObject(object)
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
    object <- minfi:::.getAnnotationObject(object)
    available <- availableAnnotation(object)
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
        get(wh, envir = object@data)
    }))
    if(!is.null(lociNames))
        out <- out[rownames(out) %in% lociNames,]
    if(dropNonMapping) {
        seqOrder <- minfi:::.seqnames.order
        wh <- which(! out$chr %in% seqOrder)
        out <- out[wh,]
    } else {
        seqOrder <- minfi:::.seqnames.order.all
    }
    if(orderByLocation) {
        sp <- split(out, out$chr)
        sp <- sp[seqOrder[seqOrder %in% names(sp)]]
        out <- do.call(rbind, lapply(sp, function(df) {
            od <- order(df$pos)
            df[od,]
        }))
    }
    out
}
