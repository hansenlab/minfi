# Author: Jean-Philippe Fortin
# May 6th 2015

# Example

## library(minfiData)
## library(lineprof)
## GMset <- mapToGenome(MsetEx)
## prof1 <- lineprof({
## con <- createCorMatrix(GMset, res=500*1000)
## })
## prof2 <- lineprof({
## ab <- extractAB(con)
## })
## prof3 <- lineprof({
## a <- compartments(GMset, resolution=500*1000)
## })

compartments <- function(object, resolution = 100*1000, what="OpenSea",
                         chr = "chr22", method = c("pearson", "spearman"),
                         keep = TRUE){
    .isGenomic(object)
    stopifnot(length(chr) == 1 && chr %in% seqlevels(object))
    method <- match.arg(method)
    gr <- createCorMatrix(object = object, resolution = resolution,
                           what = what, chr = chr, method = method)
    gr <- extractAB(gr, keep = keep)
    gr$compartment <- .extractOpenClosed(gr)
    gr
}

createCorMatrix <- function(object, resolution = 100*1000, what = "OpenSea",
                             chr = "chr22", method = c("pearson", "spearman")) {
    .isGenomic(object)
    stopifnot(length(chr) == 1 && chr %in% seqlevels(object))
    method <- match.arg(method)

    if(is(object, "GenomicMethylSet"))
        object <- ratioConvert(object, what = "M", keepCN = FALSE)
    assay(object, "M") <- .imputeMatrix(getM(object))
    ## Next we subset to a chromosome, keep OpenSea probes and remove SNPs
    seqlevels(object, force = TRUE) <- chr
    object <- object[getIslandStatus(object) %in% what,]
    object <- dropLociWithSnps(object, snps = c("CpG", "SBE"), maf = 0.01)

    gr.unbinnedCor <- granges(object)
    gr.unbinnedCor$cor.matrix <- cor(t(getM(object)), method = method)
    gr.cor <- .returnBinnedMatrix(gr.unbinnedCor, resolution = resolution)
    gr.cor <- .removeBadBins(gr.cor)
    gr.cor
}

.imputeMatrix <- function(matrix){
    matrix[is.infinite(matrix) & matrix >0] <- max(matrix[is.finite(matrix)])
    matrix[is.infinite(matrix) & matrix <0] <- min(matrix[is.finite(matrix)])
    
    ## Imputation of the missing values:
    missing <- which(is.na(matrix), arr.ind=TRUE)
    if (length(missing)!=0){
        for (j in 1:nrow(missing)){
            mean <- mean(matrix[missing[j,1],][is.finite(matrix[missing[j,1],])], na.rm=TRUE)
            matrix[missing[j,1],missing[j,2]] <- mean
        }
    }
    matrix
}

.returnBinnedMatrix <- function(gr.unbinnedCor, resolution){
    
    bin2D <- function(matrix, ids, n){
        unique.ids <- sort(unique(ids))
        bin.matrix <- matrix(0,n,n)
        m <- length(unique.ids)
        for (i in 1:n){
            for (j in 1:n){
                ## We should be able to speed this one up a lot
                indices1 <- which(ids==unique.ids[i])
                indices2 <- which(ids==unique.ids[j])
                bin.matrix[unique.ids[i],unique.ids[j]] <- median(matrix[indices1,indices2],na.rm=TRUE)
            }
        }
        bin.matrix[is.na(bin.matrix)] <- 1 # Should not be necessary... 
        return(bin.matrix)
    }

    ## FIXME
    chr.lengths <- structure(c(249250621L, 243199373L, 198022430L, 191154276L, 180915260L, 
                               171115067L, 159138663L, 146364022L, 141213431L, 135534747L, 135006516L, 
                               133851895L, 115169878L, 107349540L, 102531392L, 90354753L, 81195210L, 
                               78077248L, 59128983L, 63025520L, 48129895L, 51304566L, 155270560L, 
                               59373566L),
                             .Names = c("chr1", "chr2", "chr3", "chr4", "chr5", 
                             "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                             "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
                             "chr21", "chr22", "chrX", "chrY"))
    seqlengths(gr.unbinnedCor) <- chr.lengths[seqlevels(gr.unbinnedCor)]

    stopifnot(length(seqlevels(gr.unbinnedCor)) == 1 && !is.na(seqlengths(gr.unbinnedCor)))
    stopifnot("cor.matrix" %in% names(mcols(gr.unbinnedCor)))

    gr.binnedCor <- tileGenome(seqlengths = seqlengths(gr.unbinnedCor),
                              tilewidth = resolution, cut.last.tile.in.chrom = TRUE)
    ids <- subjectHits(findOverlaps(gr.unbinnedCor, gr.binnedCor))
    gr.binnedCor$cor.matrix <- bin2D(gr.unbinnedCor$cor.matrix, ids, n = length(gr.binnedCor))
    gr.binnedCor
}

.removeBadBins <- function(gr){
    n <- nrow(gr$cor.matrix)
    good.bins  <- which(colSums(gr$cor.matrix==0) != n)
    if(length(good.bins) < n) {
        gr <- gr[good.bins]
        gr$cor.matrix <- gr$cor.matrix[, good.bins]
    }
    return(gr)
}

extractAB <- function(gr, keep = TRUE, svdMethod = "qr"){
    if (! (is(gr, "GRanges") && "cor.matrix" %in% names(mcols(gr)))) {
        stop("'gr' must be an object created by createCorMatrix")
    }
    pc <- .getFirstPC(gr$cor.matrix, method = svdMethod)
    pc <- .meanSmoother(pc)
    pc <- .unitarize(pc)
    ## Fixing sign of eigenvector
    if (cor(colSums(gr$cor.matrix), pc) <0 ){
        pc <- -pc
    }
    pc <- pc*sqrt(length(pc))
    gr$pc <- pc
    
    if (!keep) {
        gr$cor.matrix <- NULL
    }
    return(gr)
}

.extractOpenClosed <- function(gr, cutoff = 0){
    pc <- gr$pc
    ifelse(pc < cutoff, "open", "closed")
}

.getFirstPC <- function(matrix, method){
    ## Centering the matrix
    center <- rowMeans(matrix, na.rm = TRUE)
    matrix <- sweep(matrix, 1L, center, check.margin = FALSE)
    ## if(method == "nipals")
    ##     pc <- mixOmics::nipals(matrix, ncomp = 1)$p[,1]
    pc <- .fsvd(matrix, k = 1, method = method)$u
    return(pc)
}


.meanSmoother <- function(x, k=1, iter=2, na.rm=TRUE){
    meanSmoother.internal <- function(x, k=1, na.rm=TRUE){
        n <- length(x)
        y <- rep(NA,n)
        
        window.mean <- function(x, j, k, na.rm=na.rm){
            if (k>=1){
                return(mean(x[(j-(k+1)):(j+k)], na.rm=na.rm))
            } else {
                return(x[j])
            }    
        }
        
        for (i in (k+1):(n-k)){
            y[i] <- window.mean(x,i,k, na.rm)
        }
        for (i in 1:k){
            y[i] <- window.mean(x,i,i-1, na.rm)
        }
        for (i in (n-k+1):n){
            y[i] <- window.mean(x,i,n-i,na.rm)
        }
        y
    }
    
    for (i in 1:iter){
        x <- meanSmoother.internal(x, k=k, na.rm=na.rm)
    }
    x
}

.unitarize <- function (x, medianCenter = TRUE) {
    if(medianCenter) {
        x <- x - median(x, na.rm = TRUE)
    }
    bad <- is.na(x)
    x[!bad] <- x[!bad] / sqrt(sum(x[!bad]^2))
    n.bad <- sum(bad)
    if (n.bad > 0){
        cat(sprintf("[.unitarize] %i missing values were ignored.\n", n.bad))
    }
    return(x)
}

.fsvd <- function(A, k, i = 1, p = 2, method = c("qr", "svd", "exact")){
    method <- match.arg(method)
    l <- k + p 
    n <- ncol(A)
    m <- nrow(A)
    tall <- TRUE
    ## In the case a WIDE matrix is provided:
    if (m < n){
        A <- t(A) # This is potentially expensive
        n <- ncol(A)
        m <- nrow(A)
        tall <- FALSE
    }
    if (l > ncol(A) && method != "exact"){
        stop("(k+p) is greater than the number of columns. Please decrease the value of k.")
    }
    ## Construct G to be n x l and Gaussian
    G <- matrix(rnorm(n*l,0,1), nrow=n, ncol=l)
    ## SVD approach
    if (method == "svd"){
        ## Power method: 
        H <- A %*% G # m x l matrix
        for (j in 1:i){
            H <- A %*% (crossprod(A, H))
        }
        ## We use a SVD to find an othogonal basis Q:
        ## H = F %*% Omega %*% t(S)
        svd <- svd(crossprod(H))
        F   <- svd$u # l x l
        omega <- diag(1/sqrt(svd$d)) # l x l
        S <- H %*% F %*% omega # m x l 
        ## Define the orthogonal basis:
        Q <- S[,1:k,drop=FALSE] # m x k
        ## T <- t(A) %*% Q # n x k 
        ## T <- t(T)
        T <- crossprod(Q, A)
    } 
    ## QR approach
    if (method == "qr"){
        ## Need to create a list of H matrices
        h.list <- vector("list", i+1)
        h.list[[1]] <- A %*% G 
        for (j in 2:(i+1)){
            h.list[[j]] <- A %*% (crossprod(A, h.list[[j-1]]))
        }
        H <- do.call("cbind",h.list) # n x [(1+1)l] matrix
        ## QR algorithm
        Q <- qr.Q(qr(H,0))
        ## T <- t(A)%*%Q # n x [(i+1)l]
        ## T <- t(T)
        T <- crossprod(Q, A)
    }
    if (method == "svd" | method == "qr"){
        svd <- svd(T)
        u <- Q %*% svd$u[,1:k,drop=FALSE]
        v <- svd$v[,1:k,drop=FALSE]
        d <- svd$d[1:k]
    } else { # For exact SVD:
        svd <- svd(A)
        u <- svd$u[,1:k,drop=FALSE]
        v <- svd$v[,1:k,drop=FALSE]
        d <- svd$d[1:k]
    }
    if (!tall){
        uu <- v
        v <- u
        u <- uu
    } 
    list(u=u, v=v, d=d)
}
