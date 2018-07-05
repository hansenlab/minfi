# Internal functions -----------------------------------------------------------

.imputeMatrix <- function(matrix, lower = min(matrix[is.finite(matrix)]),
                          upper = max(matrix[is.finite(matrix)])) {
    # NOTE: Does not support DelayedMatrix
    stopifnot(is.matrix(matrix))

    # Replace positive (resp. negative) infinite values by `upper` (resp.
    # `lower`).
    matrix <- pmax(pmin(matrix, upper), lower)

    # Imputation of the missing values: Replace NA/NaN by (finite) row mean
    rows_with_NAs <- rowAnyNAs(matrix)
    matrix[rows_with_NAs, ] <- rowMeans2(
        x = matrix,
        rows = rows_with_NAs,
        na.rm = TRUE)

    matrix
}

.returnBinnedMatrix <- function(gr.unbinnedCor, resolution) {

    bin2D <- function(matrix, ids, n) {
        unique.ids <- sort(unique(ids))
        bin.matrix <- matrix(0, n, n)
        m <- length(unique.ids)
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                # TODO: We should be able to speed this one up a lot
                indices1 <- which(ids == unique.ids[i])
                indices2 <- which(ids == unique.ids[j])
                bin.matrix[unique.ids[i], unique.ids[j]] <- median(
                    x = matrix[indices1, indices2],
                    na.rm = TRUE)
            }
        }
        # TODO: This line should not be necessary
        bin.matrix[is.na(bin.matrix)] <- 1
        bin.matrix
    }

    chr.lengths <- structure(
        c(249250621L, 243199373L, 198022430L, 191154276L, 180915260L,
          171115067L, 159138663L, 146364022L, 141213431L, 135534747L,
          135006516L, 133851895L, 115169878L, 107349540L, 102531392L,
          90354753L, 81195210L, 78077248L, 59128983L, 63025520L, 48129895L,
          51304566L, 155270560L, 59373566L),
        .Names = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                   "chr21", "chr22", "chrX", "chrY"))
    seqlengths(gr.unbinnedCor) <- chr.lengths[seqlevels(gr.unbinnedCor)]

    stopifnot(length(seqlevels(gr.unbinnedCor)) == 1 &&
                  !is.na(seqlengths(gr.unbinnedCor)))
    stopifnot("cor.matrix" %in% names(mcols(gr.unbinnedCor)))
    stopifnot(is.matrix(gr.unbinnedCor$cor.matrix))

    gr.binnedCor <- tileGenome(
        seqlengths = seqlengths(gr.unbinnedCor),
        tilewidth = resolution,
        cut.last.tile.in.chrom = TRUE)
    ids <- subjectHits(findOverlaps(gr.unbinnedCor, gr.binnedCor))
    gr.binnedCor$cor.matrix <- bin2D(
        matrix = gr.unbinnedCor$cor.matrix,
        ids = ids,
        n = length(gr.binnedCor))
    gr.binnedCor
}

.removeBadBins <- function(gr) {
    n <- length(gr)
    good.bins <- !colAlls(gr$cor.matrix, value = 0)
    gr <- gr[good.bins]
    gr$cor.matrix <- gr$cor.matrix[, good.bins]
    gr
}

.extractOpenClosed <- function(gr, cutoff = 0){
    pc <- gr$pc
    ifelse(pc < cutoff, "open", "closed")
}

.getFirstPC <- function(matrix, method) {
    stopifnot(is.matrix(matrix))

    # Centre the rows of the matrix
    center <- rowMeans2(matrix, na.rm = TRUE)
    matrix <- sweep(
        x = matrix,
        MARGIN = 1L,
        STATS = center,
        FUN = "-",
        check.margin = FALSE)
    .fsvd(matrix, k = 1, method = method)$u
}

.meanSmoother <- function(x, k = 1L, iter = 2L, na.rm = TRUE) {
    meanSmoother.internal <- function(x, k = 1L, na.rm = TRUE) {
        n <- length(x)
        y <- rep(NA_real_, n)

        window.mean <- function(x, j, k, na.rm = na.rm){
            if (k >= 1) {
                return(mean(x[seq(j - (k + 1L), j + k)], na.rm = na.rm))
            } else {
                x[j]
            }
        }

        for (i in (seq(k + 1L, n - k))) {
            y[i] <- window.mean(x, i, k, na.rm)
        }
        for (i in seq_len(k)) {
            y[i] <- window.mean(x, i, i - 1L, na.rm)
        }
        for (i in seq(n - k + 1L, n)) {
            y[i] <- window.mean(x, i, n - i, na.rm)
        }
        y
    }

    for (i in seq_len(iter)) {
        x <- meanSmoother.internal(x, k = k, na.rm = na.rm)
    }
    x
}

.unitarize <- function(x, medianCenter = TRUE) {
    if (medianCenter) x <- x - median(x, na.rm = TRUE)
    bad <- is.na(x)
    x[!bad] <- x[!bad] / sqrt(sum(x[!bad]^2))
    n.bad <- sum(bad)
    if (n.bad > 0) {
        message(
            sprintf("[.unitarize] %i missing values were ignored.\n", n.bad))
    }
    x
}

.fsvd <- function(A, k, i = 1, p = 2, method = c("qr", "svd", "exact")) {
    stopifnot(is.matrix(A))
    method <- match.arg(method)
    l <- k + p
    n <- ncol(A)
    m <- nrow(A)
    tall <- TRUE
    # NOTE: In the case a WIDE matrix is provided:
    if (m < n) {
        A <- t(A) # This is potentially expensive
        n <- ncol(A)
        m <- nrow(A)
        tall <- FALSE
    }
    if (l > ncol(A) && method != "exact") {
        stop("(k+p) is greater than the number of columns. Please decrease ",
             "the value of k.")
    }

    # Construct G to be n x l and Gaussian
    G <- matrix(rnorm(n * l, 0, 1), nrow = n, ncol = l)

    if (method == "svd") {
        # Power method:
        H <- A %*% G # m x l matrix
        for (j in seq_len(i)) {
            H <- A %*% (crossprod(A, H))
        }
        # NOTE: We use a SVD to find an othogonal basis Q:
        svd <- svd(crossprod(H))
        FF <- svd$u # l x l
        omega <- diag(1/sqrt(svd$d)) # l x l
        S <- H %*% FF %*% omega # m x l
        # Define the orthogonal basis:
        Q <- S[, seq_len(k), drop = FALSE] # m x k
        TT <- crossprod(Q, A)
    } else if (method == "qr") {
        # NOTE: Need to create a list of H matrices
        h.list <- vector("list", i + 1L)
        h.list[[1]] <- A %*% G
        for (j in seq(2, i + 1L)) {
            h.list[[j]] <- A %*% (crossprod(A, h.list[[j - 1L]]))
        }
        H <- do.call("cbind", h.list) # n x [(1+1)l] matrix
        # QR algorithm
        Q <- qr.Q(qr(H, 0))
        TT <- crossprod(Q, A)
    }
    if (method == "svd" || method == "qr") {
        svd <- svd(TT)
        u <- Q %*% svd$u[,seq_len(k), drop = FALSE]
        v <- svd$v[, seq_len(k), drop = FALSE]
        d <- svd$d[seq_len(k)]
    } else {
        # Exact SVD
        svd <- svd(A)
        u <- svd$u[, seq_len(k), drop = FALSE]
        v <- svd$v[, seq_len(k), drop = FALSE]
        d <- svd$d[seq_len(k)]
    }
    if (!tall) {
        uu <- v
        v <- u
        u <- uu
    }
    list(u = u, v = v, d = d)
}

# Exported functions -----------------------------------------------------------

compartments <- function(object, resolution = 100000L, what = "OpenSea",
                         chr = "chr22", method = c("pearson", "spearman"),
                         keep = TRUE) {

    # Check inputs
    .isMatrixBackedOrWarning(object, FUN = "compartments")
    .isGenomicOrStop(object)
    stopifnot(length(chr) == 1 && chr %in% seqlevels(object))
    method <- match.arg(method)

    # Construct correlation matrix
    gr <- createCorMatrix(
        object = object,
        resolution = resolution,
        what = what,
        chr = chr,
        method = method)

    # Extract A/B compartments
    gr <- extractAB(gr, keep = keep)
    gr$compartment <- .extractOpenClosed(gr)

    gr
}

createCorMatrix <- function(object, resolution = 100 * 1000, what = "OpenSea",
                            chr = "chr22", method = c("pearson", "spearman")) {
    # Check inputs
    .isMatrixBackedOrWarning(object, "createCorMatrix")
    .isGenomicOrStop(object)
    stopifnot(length(chr) == 1 && chr %in% seqlevels(object))
    method <- match.arg(method)

    if (is(object, "GenomicMethylSet")) {
        object <- ratioConvert(object, what = "M", keepCN = FALSE)
    }

    # Compute lower and upper bounds used in imputation
    # NOTE: In the original implementation these bounds are computed from the
    #       complete `M`, prior to subsetting.
    # NOTE: `assay(object, "M", withDimnames = FALSE)` is more efficient than
    #       `getM(object)` since it uses `withDimnames = FALSE`.
    M <- assay(object, "M", withDimnames = FALSE)
    # TODO: Would be nice to use range(M, finite = TRUE) but
    #       range,DelayedArray-method does not support the `finite` argument
    #       (https://github.com/Bioconductor/DelayedArray/issues/18).
    if (is.matrix(M)) {
        M_finite_range <- range(M, finite = TRUE)
    } else if (is(M, "DelayedMatrix")) {
        M_finite_range <- range(
            unlist(
                blockApply(M, range, finite = TRUE, BPPARAM = SerialParam()),
                use.names = FALSE))
    }
    M_lower <- M_finite_range[1L]
    M_upper <- M_finite_range[2L]

    # Subset
    object <- keepSeqlevels(object, value = chr, pruning.mode = "coarse")
    object <- object[getIslandStatus(object) %in% what, ]
    object <- dropLociWithSnps(object, snps = c("CpG", "SBE"), maf = 0.01)

    # Realize subsetted M in-memory and impute
    M <- as.matrix(assay(object, "M", withDimnames = FALSE))
    imputed_M <- .imputeMatrix(matrix = M, lower = M_lower, upper = M_upper)

    # Construct unbinned object
    gr.unbinnedCor <- granges(object)
    gr.unbinnedCor$cor.matrix <- cor(t(imputed_M), method = method)

    # Construct binned object
    gr.cor <- .returnBinnedMatrix(gr.unbinnedCor, resolution = resolution)
    .removeBadBins(gr.cor)
}

extractAB <- function(gr, keep = TRUE, svdMethod = "qr"){
    if (!(is(gr, "GRanges") && "cor.matrix" %in% names(mcols(gr)))) {
        stop("'gr' must be an object created by createCorMatrix")
    }
    # Compute, smooth, and unitarize first principal component of correlation
    # matrix
    pc <- .getFirstPC(gr$cor.matrix, method = svdMethod)
    pc <- .meanSmoother(pc)
    pc <- .unitarize(pc)

    # Standardize first principal component
    if (cor(colSums(gr$cor.matrix), pc) < 0 ) {
        pc <- -pc
    }
    pc <- pc * sqrt(length(pc))
    gr$pc <- pc

    if (!keep) gr$cor.matrix <- NULL
    gr
}
