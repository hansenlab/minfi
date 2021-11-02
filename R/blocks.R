# Internal functions -----------------------------------------------------------

clusterMaker4Blocks <- function(gr, relationToIsland, islandName, maxGap,
                                maxClusterWidth) {
    # Get middle position of each island
    if (!is.character(relationToIsland) ||
        !all(unique(relationToIsland) %in%
             c("OpenSea", "Island", "Shelf", "N_Shelf", "S_Shelf", "Shore",
               "N_Shore", "S_Shore")) ||
        length(gr) != length(relationToIsland)) {
        stop("argument 'relationToIsland' is either not a character or seems ",
             "to have wrong values")
    }
    if (!is.character(islandName) || length(gr) != length(islandName)) {
        stop("argument 'islandName' is not a character or it has the wrong ",
             "length")
    }

    fullName <- paste0(islandName, relationToIsland)
    isNotSea <- (fullName != "OpenSea")

    pos <- start(gr)
    # NOTE: These are the "average positions" on shelfs, shores, and islands
    ipos <- round(tapply(pos[isNotSea], fullName[isNotSea], mean))
    # Make non-island/shore/shelf positions their own island
    islFactor <- ifelse(isNotSea, fullName, seq_along(pos))
    # Check which of the new positions correspond to islands/shores/shelfs,
    # i.e. open sea. Assign probes in islands/shores/shelves just one position
    pos[isNotSea] <- ipos[islFactor[isNotSea]]
    # make first pass clusters with new positions
    pns <- boundedClusterMaker(
        chr = as.numeric(seqnames(gr)),
        pos = pos,
        assumeSorted = TRUE,
        maxGap = maxGap,
        maxClusterWidth = maxClusterWidth)
    # But islands must be on their own so we change those beginning and ends of
    # each genomic region type: island, shore, shelf
    types <- unique(relationToIsland)
    # Take out open sea
    types <- types[types != "OpenSea"]
    for (i in types) {
        ind <- relationToIsland == i
        StartEnd <- abs(diff(c(0,ind) != 0))
        add2pns <- cumsum(StartEnd)
        pns <- pns + add2pns
    }
    # where are the islands?
    pns <- as.numeric(as.factor(pns))

    # annotation
    type <- sub("^[NS]_", "", relationToIsland)

    data.frame(pns = pns, type = type)
}


cpgCollapseAnnotation <- function(gr, relationToIsland, islandName,
                                  maxGap = 500, maxClusterWidth = 1500,
                                  blockMaxGap = 2.5*10^5, verbose = TRUE) {

    # Check inputs
    if (is.unsorted(order(gr))) stop("object has to be ordered.")
    # NOTE: This function aggregates probes in islands and the rest are
    #       aggregated into clusters. Then annotation is created for the new
    #       aggregated data the indexes of the original probes along with the
    #       group ids are returned indexes and pns respectively. Note that
    #       these two are redudant but they save me work
    # Make sure chr is a factor.
    # TODO (Kasper): Should this be done from get-go?

    if (verbose) {
        message("[cpgCollapseAnnotation] Clustering islands and clusters of ",
                "probes.\n")
    }

    # Block tab has the block group ids
    blocktab <- clusterMaker4Blocks(
        gr = gr,
        relationToIsland = relationToIsland,
        islandName = islandName,
        maxGap = maxGap,
        maxClusterWidth = maxClusterWidth)

    # Split rows by group id
    groupIndexes <- split(seq_along(blocktab$pns), blocktab$pns)

    # make an object with results
    if (verbose) message("[cpgCollapseAnnotation] Computing new annotation.\n")
    tmpRanges <- t(sapply(groupIndexes, function(ind) range(start(gr)[ind])))

    anno <- GRanges(
        seqnames = Rle(
            tapply(as.vector(seqnames(gr)), blocktab$pns, function(x) x[1])),
        ranges = IRanges(start = tmpRanges[,1], end = tmpRanges[,2]),
        id = as.numeric(names(groupIndexes)),
        type = as.vector(
            tapply(
                as.character(blocktab$type), blocktab$pns, function(x) x[1])))
    res <- list(anno = anno, indexes = groupIndexes)
    seql <- seqlevels(res$anno)
    seqlevels(res$anno, pruning.mode = "coarse") <-
        .seqnames.order[.seqnames.order %in% seql]

    if (verbose) message("[cpgCollapseAnnotation] Defining blocks.\n")
    ind <- (res$anno$type == "OpenSea")
    pns <- rep(NA, length(res$anno))
    pns[ind] <- clusterMaker(
        chr = as.numeric(seqnames(res$anno[ind,])),
        pos = start(res$anno[ind,]),
        maxGap = blockMaxGap)
    res$anno$blockgroup <- pns
    res$pns <- blocktab$pns
    res
}

# Internal generics ------------------------------------------------------------

# `x` is either `meth_signal` or `cn`
# `...` are additional arguments passed to methods.
setGeneric(
    ".cpgCollapse",
    function(x, Indexes, dataSummary, na.rm, verbose, ...)
        standardGeneric(".cpgCollapse"),
    signature = "x")

# Internal methods -------------------------------------------------------------

setMethod(
    ".cpgCollapse",
    "matrix",
    function(x, Indexes, dataSummary, na.rm, verbose) {

        # Set up return matrix
        n_clusters <- length(Indexes)
        n_cpgs_per_cluster <- lengths(Indexes)
        stopifnot(all(n_cpgs_per_cluster) > 0)
        # NOTE: .cpgCollapse() can return non-integer values, so fill a numeric
        #        matrix
        collapsed_x <- matrix(NA_real_, nrow = n_clusters, ncol = ncol(x))

        # Process clusters with a single CpG
        cluster_idx <- n_cpgs_per_cluster == 1L
        cpg_idx <- unlist(Indexes[cluster_idx], use.names = FALSE)
        collapsed_x[cluster_idx, ] <- x[cpg_idx, , drop = FALSE]

        # Process clusters with 0 or > 1 CpG
        cluster_idx <- which(n_cpgs_per_cluster != 1L)
        for (i in cluster_idx) {
            if (verbose) if (runif(1) < 0.0001) cat(".")
            # TODO: If willing to assume dataSummary comes from
            #       DelayedMatrixStats, could use `rows` rather than explicit
            #       subsetting.
            collapsed_x[i, ] <- dataSummary(
                x[Indexes[[i]], , drop = FALSE],
                na.rm = na.rm)
        }
        if (verbose) cat("\n")
        collapsed_x
    })

setMethod(
    ".cpgCollapse",
    "DelayedMatrix",
    function(x, Indexes, dataSummary, na.rm, verbose, BPREDO = list(),
             BPPARAM = SerialParam()) {

        # Set up intermediate RealizationSink object of appropriate dimensions
        # and type
        n_clusters <- length(Indexes)
        n_cpgs_per_cluster <- lengths(Indexes)
        # NOTE: This is ultimately coerced to the output DelayedMatrix object
        # NOTE: .cpgCollapse() can return non-integer values, so fill a "double"
        #       sink
        ans_type <- "double"
        sink <- DelayedArray::AutoRealizationSink(
            dim = c(n_clusters, ncol(x)),
            type = ans_type)
        on.exit(close(sink))

        # Set up ArrayGrid instances over `x` as well as "parallel" ArrayGrid
        # instance over `sink`
        x_grid <- colAutoGrid(x)
        sink_grid <- RegularArrayGrid(
            refdim = dim(sink),
            spacings = c(nrow(sink), ncol(sink) / length(x_grid)))
        # Sanity check ArrayGrid objects have the same dim
        stopifnot(dim(x_grid) == dim(sink_grid))

        # Loop over column-blocks of `x`, perform SWAN normalization, and write
        # to `normalized_x_sink`
        blockApplyWithRealization(
            x = x,
            FUN = .cpgCollapse,
            Indexes = Indexes,
            dataSummary = dataSummary,
            na.rm = na.rm,
            verbose = verbose,
            sink = sink,
            x_grid = x_grid,
            sink_grid = sink_grid,
            BPREDO = BPREDO,
            BPPARAM = BPPARAM)

        # Return as DelayedMatrix
        as(sink, "DelayedArray")
    }
)

# Exported functions -----------------------------------------------------------

# NOTE: Collapses a minfi object into islands, shores, and defines block regions
cpgCollapse <- function(object, what = c("Beta", "M"), maxGap = 500,
                        blockMaxGap = 2.5 * 10^5, maxClusterWidth = 1500,
                        dataSummary = colMeans, na.rm = FALSE,
                        returnBlockInfo = TRUE, islandAnno = NULL,
                        verbose = TRUE, ...) {

    # Check inputs
    # TODO: ?cpgCollapse suggests `object` needn't be a
    #       Genomic[MethlSet|RatioSet] but the code assumes `granges(object)`
    #       works
    what <- match.arg(what)

    # Construct annotation
    if (verbose) message("[cpgCollapse] Creating annotation.\n")
    islands <- .getIslandAnnotation(object = object, islandAnno = islandAnno)
    relationToIsland <- islands$Relation_to_Island
    islandName <- islands$Islands_Name
    gr <- granges(object)
    anno <- cpgCollapseAnnotation(
        gr = gr,
        relationToIsland = relationToIsland,
        islandName = islandName,
        maxGap = maxGap,
        blockMaxGap = blockMaxGap,
        maxClusterWidth = maxClusterWidth,
        verbose = verbose)
    Indexes <- split(seq_along(anno$pns), anno$pns)

    # Collapse data
    if (verbose) message("[cpgCollapse] Collapsing data")
    meth_signal <- getMethSignal(object, what = what, ...)
    collapsed_meth_signal <- .cpgCollapse(
        x = meth_signal,
        Indexes = Indexes,
        dataSummary = dataSummary,
        na.rm = na.rm,
        verbose = verbose)
    cn <- getCN(object,...)
    if (!is.null(cn)) {
        collapsed_cn <- .cpgCollapse(
            x = cn,
            Indexes = Indexes,
            dataSummary = dataSummary,
            na.rm = na.rm,
            verbose = verbose)
    }

    # Construct output
    preproc <- c(collapse = "cpgCollapse", preprocessMethod(object))
    if (what == "M") {
        ret <- GenomicRatioSet(
            gr = anno$anno,
            Beta = NULL,
            M = collapsed_meth_signal,
            CN = collapsed_cn,
            colData = colData(object),
            annotation = annotation(object),
            preprocessMethod = preproc)
    }
    else {
        ret <- GenomicRatioSet(
            gr = anno$anno,
            Beta = collapsed_meth_signal,
            M = NULL,
            CN = collapsed_cn,
            colData = colData(object),
            annotation = annotation(object),
            preprocessMethod = preproc)
    }
    # NOTE: Take out annotation as we already kept it
    anno <- anno[2:3]
    if (returnBlockInfo) {
        return(list(object = ret, blockInfo = anno))
    }
    ret
}

# NOTE: blockFinder() just uses a cluster object or granges()$blockgroup
#       where clusters are constructed by
#           cpgCollpase() which calls
#               cpgCollapseAnnotation() which calls
#                   clusterMaker4Blocks() which calls
#                       boundedClusterMaker()
blockFinder <- function(object, design, coef = 2, what = c("Beta", "M"),
                        cluster = NULL, cutoff = NULL,
                        pickCutoff = FALSE, pickCutoffQ = 0.99,
                        nullMethod = c("permutation","bootstrap"),
                        smooth = TRUE, smoothFunction = locfitByCluster,
                        B = ncol(permutations), permutations = NULL,
                        verbose = TRUE, bpSpan = 2.5*10^5, ...) {

    # Check inputs
    .isMatrixBackedOrStop(object)
    if (!is(object,"GenomicRatioSet")) stop("object must be 'GenomicRatioSet'")
    if (is.null(cluster)) cluster <- granges(object)$blockgroup
    if (is.null(cluster)) stop("need 'cluster'")
    what <- match.arg(what)
    nullMethod <- match.arg(nullMethod)
    idx <- which(granges(object)$type == "OpenSea")
    if (length(idx) == 0) stop("need OpenSea types in granges(object)")


    pos <- start(granges(object)) / 2 + end(granges(object)) / 2
    res <- bumphunterEngine(
        mat = getMethSignal(object, what)[idx,],
        design = design,
        coef = coef,
        chr = as.character(seqnames(object))[idx],
        pos = pos[idx], cluster = cluster[idx],
        cutoff = cutoff,  pickCutoff = pickCutoff,
        pickCutoffQ = pickCutoffQ,
        nullMethod = nullMethod,
        smooth = smooth,
        smoothFunction = smoothFunction,
        B = B,
        permutations = permutations,
        verbose = verbose,
        bpSpan = bpSpan,...)

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

