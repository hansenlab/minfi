# AUthor: Jean-Philippe Fortin
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
    method <- match.arg(method)

    if(is(object, "GenomicMethylSet"))
        object <- ratioConvert(object, what = "M", keepCN = FALSE)
    assay(object, "M") <- .imputeMatrix(getM(object))
    ## Next we subset to a chromosome, keep OpenSea probes and remove SNPs
    seqlevels(object, force = TRUE) <- chr
    object <- object[getIslandStatus(object) %in% what,]
    object <- dropLociWithSnps(object, snps = c("CpG", "SBE"), maf = 0.01)
    
    matrix <- getM(object)
    ann <- data.frame(chr=seqnames(object), pos=start(object))
    rownames(ann) <- rownames(matrix)
    unbinnedCor <- cor(t(matrix), method = method)
    gr.cor <- .returnBinnedMatrix(matrix=unbinnedCor, ann=ann, res=resolution)
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

.removeSnps <- function(matrix){
    ## Hardcoding annotation ... also do we have existing functions for this?
    ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    ann <- ann[match(rownames(matrix), rownames(ann)),]
    snp.info <- ann[, c("CpG_maf","SBE_maf" )]
    indices1 <- which(snp.info[,"CpG_maf"] >0.01)
    indices2 <- which(snp.info[,"SBE_maf"] >0.01)
    indices3 <- which(grepl("rs", rownames(matrix)))
    snps <- union(indices1, indices2)
    if (length(snps)>0){
        matrix <- matrix[-snps,]
    }
    matrix
}

.returnBinnedMatrix <- function(matrix, ann, res){
    ## for convenience.. (UCSC hg19)
    ## FIMXE: this is ugly; can we do this other way?
    chr.lengths <- structure(c(249250621L, 243199373L, 198022430L, 191154276L, 180915260L, 
                               171115067L, 159138663L, 146364022L, 141213431L, 135534747L, 135006516L, 
                               133851895L, 115169878L, 107349540L, 102531392L, 90354753L, 81195210L, 
                               78077248L, 59128983L, 63025520L, 48129895L, 51304566L, 155270560L, 
                               59373566L),
                             .Names = c("chr1", "chr2", "chr3", "chr4", "chr5", 
                             "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                             "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
                             "chr21", "chr22", "chrX", "chrY"))
    
    bin2D <- function(matrix,ids,n){
        unique.ids <- sort(unique(ids))
        new.matrix <- matrix(0,n,n)
        m <- length(unique.ids)
        for (i in 1:n){
            for (j in 1:n){
                ## We should be able to speed this one up a lot
                indices1 <- which(ids==unique.ids[i])
                indices2 <- which(ids==unique.ids[j])
                new.matrix[unique.ids[i],unique.ids[j]] <- median(matrix[indices1,indices2],na.rm=TRUE)
            }
        }
        new.matrix[is.na(new.matrix)] <- 1 # Should not be necessary... 
        return(new.matrix)
    }
    
    loci <- rownames(matrix)
    ann <- ann[match(loci, rownames(ann)),] # Should be necessary
    chr <- as.character(ann$chr[1])
    chr.max <- chr.lengths[chr]
    chr.min <- 0
    start <- seq(chr.min,chr.max,res)
    end <- c(start[-1],chr.max) -1L
    
    bin.granges <- GRanges(seqnames = ann$chr[1],
                           ranges = IRanges(start = start, end = end))
    chr.granges <- GRanges(seqnames = ann$chr[1],
                           ranges = IRanges(start=ann$pos, end=ann$pos))
    ids <- subjectHits(findOverlaps(chr.granges, bin.granges))
    binned.matrix <- bin2D(matrix, ids, n=length(bin.granges))

    gr <- bin.granges
    gr$cor.matrix <- binned.matrix
    gr
}

## .getIslandStatus <- function(probes) {
##     ## FIXME: hardcoding annotation; also this function more or less exisit in minfi
##     ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
##     ann <- ann[match(probes, rownames(ann)), ]
##     ann$Relation_to_Island
## }


.removeBadBins <- function(gr){
    n <- nrow(gr$cor.matrix)
    good.bins  <- which(colSums(gr$cor.matrix==0)!=n)
    if(length(good.bins) < n) {
        gr <- gr[good.bins]
        gr$cor.matrix <- gr$cor.matrix[, good.bins]
    }
    return(gr)
}

extractAB <- function(gr, keep = TRUE){
    if (! (is(gr, "GRanges") && "cor.matrix" %in% names(mcols(gr)))) {
        stop("'gr' must be an object created by createCorMatrix")
    }
    
    pc <- .getFirstPC(gr$cor.matrix)
    pc <- .meanSmoother(pc)
    pc <- .unitarize(pc)
    ## Fixing sign of eigenvector
    if (cor(colSums(gr$cor.matrix), pc) <0 ){
        pc <- -pc
    }
    gr$pc <- pc
    if (!keep) {
        gr$cor.matrix <- NULL
    }
    return(gr)
}

.extractOpenClosed <- function(gr){
    pc <- gr$pc
    ifelse(pc < 0, "open", "closed")
}

.getFirstPC <- function(matrix, ncomp=1){
    ## Centering the matrix
    center <- rowMeans(matrix, na.rm = TRUE)
    matrix <- sweep(matrix, 1L, center, check.margin = FALSE)
    if (ncomp > 1){
        ## FIXME: shouldn't $p be subsetted with ncomp? So it just extracts a single pc
        return(mixOmics::nipals(matrix, ncomp = ncomp)$p)
    } else {
        return(mixOmics::nipals(matrix, ncomp = ncomp)$p[,1])
  }
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












