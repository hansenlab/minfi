.fixMethOutliers <- function(mat, K=-3, verbose = FALSE) {
    logMat <- log2(mat + 0.5)
    mu <- matrixStats::colMedians(logMat)
    sd <- matrixStats::colMads(logMat)
    cutoff <- K*sd + mu
    nFixes <- integer(ncol(mat))
    for(jj in 1:ncol(mat)) {
        wh <- which(logMat[, jj] < cutoff[jj])
        mat[wh, jj] <- 2^cutoff[jj]
        nFixes[jj] <- length(wh)
        if(verbose)
            cat(sprintf("[.fixMethOutliers] for sample %s, fixing %s outliers with 2^cutoff=%s\n",
                        jj, nFixes[jj], round(2^cutoff[jj],0)))
    }
    return(list(mat = mat, cutoff = cutoff, nFixes = nFixes))
}

fixMethOutliers <- function(object, K=-3, verbose = FALSE){
    .isMethyl(object)
    if(verbose) cat("[fixMethOutliers] fixing Meth channel\n")
    if(is(object, "GenomicMethylSet"))
        assay(object, "Meth") <- .fixMethOutliers(getMeth(object), K=K, verbose = verbose)$mat
    else
        assayDataElement(object, "Meth") <- .fixMethOutliers(getMeth(object), K=K, verbose = verbose)$mat
    if(verbose) cat("[fixMethOutliers] fixing Unmeth channel\n")
    if(is(object, "GenomicMethylSet"))
        assay(object, "Unmeth") <- .fixMethOutliers(getUnmeth(object), K=K, verbose = verbose)$mat
    else
        assayDataElement(object, "Unmeth") <- .fixMethOutliers(getUnmeth(object), K=K, verbose = verbose)$mat
      return(object)
}

addQC <- function(object, qc) {
    .pDataAdd(object, qc)
}

plotQC <- function(qc, badSampleCutoff = 10.5) {
    meds <- (qc$mMed + qc$uMed)/2
    whichBad <- which((meds < badSampleCutoff))
    plot(qc$mMed, qc$uMed,
         xlim = c(8,14), ylim = c(8,14), xaxt = "n", yaxt = "n",
         xlab = "Meth median intensity (log2)",
         ylab = "Unmeth median intensity (log2)",
         col = ifelse(1:nrow(qc) %in% whichBad, "red", "black"))
    axis(side = 1, at = c(9,11,13))
    axis(side = 2, at = c(9,11,13))
    ## abline(h = badSampleCutoff, lty = 2)
    ## abline(v = badSampleCutoff, lty = 2)
    abline(badSampleCutoff *2 , -1, lty = 2)
    if(length(whichBad) > 0)
        text(qc$mMed[whichBad], qc$uMed[whichBad] - 0.25,
             labels = whichBad, col = "red")
    legend("topleft", legend = c("good", "bad, with sample index"), pch = 1,
           col = c("black", "red"), bty = "n")
    invisible(NULL)
}

getQC <- function(object) {
    .isMethyl(object)
    U.medians <- log2(matrixStats::colMedians(getUnmeth(object)))
    M.medians <- log2(matrixStats::colMedians(getMeth(object)))
    df <- DataFrame(mMed = M.medians, uMed = U.medians)
    rownames(df) <- sampleNames(object)
    df
}

minfiQC <- function(object, fixOutliers=TRUE, verbose = FALSE){
    .isMethyl(object)
    subverbose <- max(as.integer(verbose) - 1L, 0L)
    if(fixOutliers){
        if(verbose) cat("[minfiQC] fixing outliers\n")
        fixMethOutliers(object, verbose = subverbose)
    }
    qc <- getQC(object)
    object <- addQC(object, qc = qc)
    if(!is(object, "GenomicMethylSet")) {
        sex <- getSex(mapToGenome(object))
    } else {
        sex <- getSex(object)
    }
    object <- addSex(object, sex = sex)
    return(list(object = object, qc = cbind(qc, sex)))
}

