# Internal functions -----------------------------------------------------------

.fixMethOutliers <- function(mat, K = -3, verbose = FALSE) {
    # Compute cutoff
    log2_mat <- log2(mat + 0.5)
    mu <- colMedians(log2_mat)
    sd <- colMads(log2_mat)
    cutoff <- 2 ^ (K * sd + mu)

    # NOTE: nFixes is only computed if verbose is TRUE because it's not
    #       otherwise needed
    if (verbose) {
        nFixes <- colSums2(sweep(mat, 2, cutoff, `<`))
        msg <- paste0(
            "[.fixMethOutliers] for sample %s, fixing %s outliers with ",
            "2^cutoff=%s\n")
        for (j in seq_len(ncol(mat))) {
            message(sprintf(msg, j, nFixes[j], round(cutoff[j], 0)))
        }
    }

    # Threshold matrix values so that all values are greater than or equal to
    # the column-specific `cutoff`
    sweep(mat, 2, cutoff, pmax2)
}

# Exported functions -----------------------------------------------------------

fixMethOutliers <- function(object, K = -3, verbose = FALSE) {
    .isMethylOrStop(object)
    if (verbose) message("[fixMethOutliers] fixing Meth channel\n")
    assay(object, "Meth") <- .fixMethOutliers(getMeth(object), K, verbose)
    if (verbose) message("[fixMethOutliers] fixing Unmeth channel\n")
    assay(object, "Unmeth") <- .fixMethOutliers(getUnmeth(object), K, verbose)
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
    abline(badSampleCutoff * 2 , -1, lty = 2)
    if (length(whichBad) > 0) {
        text(qc$mMed[whichBad], qc$uMed[whichBad] - 0.25,
             labels = whichBad, col = "red")
    }
    legend("topleft", legend = c("good", "bad, with sample index"), pch = 1,
           col = c("black", "red"), bty = "n")
    invisible(NULL)
}

getQC <- function(object) {
    .isMethylOrStop(object)
    U.medians <- log2(colMedians(getUnmeth(object)))
    M.medians <- log2(colMedians(getMeth(object)))
    df <- DataFrame(mMed = M.medians, uMed = U.medians)
    rownames(df) <- colnames(object)
    df
}

minfiQC <- function(object, fixOutliers=TRUE, verbose = FALSE){
    .isMethylOrStop(object)
    subverbose <- max(as.integer(verbose) - 1L, 0L)
    if (fixOutliers) {
        if (verbose) message("[minfiQC] fixing outliers\n")
        fixMethOutliers(object, verbose = subverbose)
    }
    qc <- getQC(object)
    object <- addQC(object, qc = qc)
    if (!is(object, "GenomicMethylSet")) {
        sex <- getSex(mapToGenome(object))
    } else {
        sex <- getSex(object)
    }
    object <- addSex(object, sex = sex)
    return(list(object = object, qc = cbind(qc, sex)))
}
