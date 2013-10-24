blockFinder <- function(object, design, coef = 2, what = c("Beta", "M"),
                        cluster = NULL, cutoff = NULL, pickCutoff = FALSE,
                        smooth = TRUE, smoothFunction = locfitByCluster,
                        B = 1000, verbose = TRUE, bpSpan = 2.5*10^5, ...) {
    if (!is(object,"GenomicRatioSet")) stop("object must be 'GenomicRatioSet'")

    if(is.null(cluster)) cluster <- granges(object)$blockgroup
    if(is.null(cluster))
        stop("need 'cluster'")

    what <- match.arg(what)
    idx <- which(granges(object)$type == "OpenSea")
    if(length(idx) == 0) stop("need OpenSea types in granges(object)")
    pos <- start(granges(object)) /2 + end(granges(object))/2
    res <- bumphunterEngine(getMethSignal(object, what)[idx,], design = design,
                            chr = as.character(seqnames(object))[idx],
                            pos = pos[idx], cluster = cluster[idx],
                            cutoff = cutoff, coef = coef, pickCutoff = pickCutoff,
                            smooth = smooth, smoothFunction = smoothFunction,
                            B = B, bpSpan = bpSpan, verbose = verbose)

    ## FIXME: reindex like below
    res$coef <- bumphunter:::.getEstimate(getMethSignal(object, what), design, coef = coef)

    ## Re-indexing because we only fit the model on the idx indexes
    res$table$indexStart <- idx[res$table$indexStart]
    res$table$indexEnd <- idx[res$table$indexEnd]
    fitted <- rep(NA, length(granges(object)))
    fitted[idx] <- res$fitted
    res$fitted <- fitted
    pvaluesMarginal <- rep(NA,length(granges(object)))
    pvaluesMarginal[idx] <- res$pvaluesMarginal
    res$pvaluesMarginal <- pvaluesMarginal

    return(res)
}



## blockFinder just uses a cluster object or granges()$blockgroup
## clusters are constructed by
##   cpgCollpase which calls
##     cpgCollapseAnnotation which calls
##       clusterMaker4Blocks which calls
##         boundedClusterMaker


##### collapses a minfi object into islands, shores, and defines block regions
cpgCollapse <- function(object, what = c("Beta", "M"), maxGap = 500,
                        blockMaxGap = 2.5*10^5, maxClusterWidth = 1500,
                        dataSummary = colMeans, na.rm = FALSE,
                        returnBlockInfo = TRUE, islandAnno = NULL, verbose = TRUE, ...) { ### ... is for illumna type
    what <- match.arg(what)
    gr <- granges(object)
    islands <- .getIslandAnnotation(object = object, islandAnno = islandAnno)
    relationToIsland <- islands$Relation_to_Island
    islandName <- islands$Islands_Name
    if(verbose) cat("[cpgCollapse] Creating annotation.\n")
    anno <- cpgCollapseAnnotation(gr, relationToIsland, islandName,
                                  maxGap = maxGap, blockMaxGap = blockMaxGap,
                                  maxClusterWidth = maxClusterWidth,
                                  verbose = verbose)
    y <- getMethSignal(object, what = what, ...)
    a <- getCN(object,...)
    
    if(verbose) cat("[cpgCollapse] Collapsing data")
    
    Indexes <- split(seq(along = anno$pns), anno$pns)
    yy  <- matrix(0, length(Indexes), ncol(y))
    aa  <- matrix(0,length(Indexes), ncol(y))
    Ns <- sapply(Indexes, length)
    ones <- (Ns == 1)
    ind <- which(ones)
    yy[ind,] <- y[unlist(Indexes[ind]),]
    aa[ind,] <- a[unlist(Indexes[ind]),]

    ind <- which(!ones)
    for(i in ind) {
        if(verbose) if(runif(1) < 0.0001) cat(".")
        yy[i,] <- dataSummary(y[Indexes[[i]],, drop=FALSE], na.rm = na.rm)
        aa[i,] <- dataSummary(a[Indexes[[i]],, drop=FALSE], na.rm = na.rm)
    }
    if(verbose) cat("\n")
    
    preproc <- c(collapse = "cpgCollapse", preprocessMethod(object))
    if(what == "M") {
        ret <- GenomicRatioSet(gr = anno$anno,
                               Beta = NULL, M = yy, CN = aa,
                               pData = pData(object),
                               annotation = annotation(object),
                               preprocessMethod = preproc) 
    }
    else {
        ret <- GenomicRatioSet(gr = anno$anno,
                               Beta = yy, M = NULL, CN = aa,
                               pData = pData(object),
                               annotation = annotation(object),
                               preprocessMethod = preproc)
    }
    anno <- anno[2:3] ## take out annotation as we already kept it
    if(returnBlockInfo) return(list(object = ret, blockInfo = anno)) else return(ret)
}

cpgCollapseAnnotation <- function(gr, relationToIsland, islandName,
                                  maxGap = 500, maxClusterWidth = 1500,
                                  blockMaxGap = 2.5*10^5, verbose = TRUE) {

    if(is.unsorted(order(gr))) stop("object has to be ordered.")
    ## this function aggregates probes in islands
    ## and the rest are aggregated into clusters
    ## then annotation is created for the new aggregated data
    ## the indexes of the original probes along with the group ids are
    ## returned indexes and pns respectively. note these two are redudant
    ## but they save me work
    
    ## make sure chr is a factor.. Kasper, should this be done from get-go?
    if(verbose) cat("[cpgCollapseAnnotation] Clustering islands and clusters of probes.\n")
    ## block tab has the block group ids
    blocktab <- clusterMaker4Blocks(gr, relationToIsland, islandName,
                                    maxClusterWidth = maxClusterWidth,
                                    maxGap = maxGap)
    
    ## split rows by group id
    groupIndexes <- split(seq(along = blocktab$pns), blocktab$pns)

    ## make an object with results
    if(verbose) cat("[cpgCollapseAnnotation] Computing new annotation.\n")
    tmpRanges <- t(sapply(groupIndexes, function(ind) range(start(gr)[ind])))

    anno <- GRanges(seqnames = Rle(tapply(as.vector(seqnames(gr)), blocktab$pns,function(x) x[1])),
                    ranges = IRanges(start = tmpRanges[,1], end = tmpRanges[,2]),
                    id = as.numeric(names(groupIndexes)),
                    type = as.vector(tapply(as.character(blocktab$type), blocktab$pns, function(x) x[1])))
    res <- list(anno = anno, indexes = groupIndexes)
    seql <- seqlevels(res$anno)
    seqlevels(res$anno, force = TRUE) <- .seqnames.order[.seqnames.order %in% seql]
    
    if(verbose) cat("[cpgCollapseAnnotation] Defining blocks.\n")
    ind <- (res$anno$type == "OpenSea")
    pns <- rep(NA, length(res$anno))
    pns[ind] <- clusterMaker(as.numeric(seqnames(res$anno[ind,])),
                             start(res$anno[ind,]), maxGap = blockMaxGap)
    res$anno$blockgroup <- pns
    res$pns <- blocktab$pns
    return(res)
}

clusterMaker4Blocks <- function(gr, relationToIsland, islandName, maxGap, maxClusterWidth) {
    ## get middle position of each island
    if(!is.character(relationToIsland) ||
       ! all(unique(relationToIsland) %in% c("OpenSea", "Island", "Shelf", "N_Shelf",
                                             "S_Shelf", "Shore", "N_Shore", "S_Shore")) ||
       length(gr) != length(relationToIsland))
        stop("argument 'relationToIsland' is either not a character or seems to have wrong values")
    if(!is.character(islandName) ||
       length(gr) != length(islandName))
        stop("argument 'islandName' is not a character or it has the wrong length")
    
    fullName <- paste0(islandName, relationToIsland)
    isNotSea <- (fullName != "OpenSea")

    pos <- start(gr)
    ## these are the "average positions" on shelfs, shores, and islands
    ipos <- round(tapply(pos[isNotSea], fullName[isNotSea], mean))
    ## make non-island/shore/shelf positions their own island
    islFactor <- ifelse(isNotSea, fullName, seq_along(pos))
    ## check which of the new positions correspond to
    ## islands/shores/shelfs, i.e. open sea
    ## assign probes in islands/shores/shelves just one position
    pos[isNotSea] <- ipos[islFactor[isNotSea]]
    ## make first pass clusters with new positions
    pns <- boundedClusterMaker(chr = as.numeric(seqnames(gr)), pos = pos,
                               assumeSorted = TRUE, maxGap = maxGap,
                               maxClusterWidth = maxClusterWidth)
    ## but islands must be on their own so we change those
    ## beginning and ends of each genomic region type:
    ## island,shore,shelf)
    types <- unique(relationToIsland)
    ## take out open sea
    types <- types[types != "OpenSea"]
    for(i in types) {
        ind <- relationToIsland == i
        StartEnd <- abs(diff(c(0,ind) != 0))
        add2pns <- cumsum(StartEnd)
        pns <- pns + add2pns
    }    
    ## where are the islands
    pns <- as.numeric(as.factor(pns))
    
    ## annotation
    type <- sub("^[NS]_", "", relationToIsland)
    return(data.frame(pns = pns, type = type))
}

