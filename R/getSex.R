getSex <- function(object = NULL, cutoff = -2){
    .isGenomic(object)
    if(is(object, "GenomicMethylSet"))
        CN <- getCN(object)
    if(is(object, "GenomicRatioSet"))
        CN <- getCN(object)
    ## FIXME: add test for logarithmic scale or non-log scale
    xIndex <- which(seqnames(object) == "chrX")
    yIndex <- which(seqnames(object) == "chrY")
    out <- .getSex(CN = CN, xIndex = xIndex,
                   yIndex = yIndex, cutoff = cutoff)
    return(out)
}

addSex <- function(object, sex = NULL) {
    if(is.null(sex))
        sex <- getSex(object)$predictedSex
    if(is(sex, "DataFrame") && "predictedSex" %in% names(sex))
        sex <- sex$predictedSex
    sex <- .checkSex(sex)
    .pDataAdd(object, DataFrame(predictedSex = sex))
}

plotSex <- function(object, id = NULL) {
    stopifnot(all(c("predictedSex", "xMed", "yMed") %in% names(object)))
    if(is.null(id))
        id <- 1:length(object$predictedSex)
    if(length(id) != length(object$predictedSex))
        stop("id length must match number of samples.")
    plot(object$xMed, object$yMed, type = "n",
         xlab = "X chr, median total intensity (log2)",
         ylab = "Y chr, median total intensity (log2)")
    text(object$xMed, object$yMed, id,
         col=ifelse(object$predictedSex == "M", "deepskyblue", "deeppink3"))
    legend("bottomleft", c("M","F"), col = c("deepskyblue", "deeppink3"), pch = 16)
}

.getSex <- function(CN = NULL, xIndex = NULL, yIndex = NULL, cutoff=-2) {
    if(is.null(CN) | is.null(xIndex) | is.null(yIndex))
        stop("must provide CN, xIndex, and yIndex")
    ## FIXME: does not handle only females or only males
    ## this ought to be handled by the 'centers' (see below) being too close together
    xMed <- matrixStats::colMedians(CN[xIndex,], na.rm=TRUE)
    yMed <- matrixStats::colMedians(CN[yIndex,], na.rm=TRUE)
    dd <- yMed - xMed
    k <- kmeans(dd, centers = c(min(dd), max(dd)))

    sex0 <- ifelse(dd < cutoff, "F", "M")
    sex0 <- .checkSex(sex0)
    sex1 <- ifelse(k$cluster == which.min(k$centers), "F", "M")
    sex1 <- .checkSex(sex1)

    if(!identical(sex0,sex1))
    warning("An inconsistency was encountered while determining sex. One possibility is that only one sex is present. We recommend further checks, for example with the plotSex function.")
    df <- DataFrame(xMed = xMed, yMed = yMed, predictedSex = sex0)
    rownames(df) <- colnames(CN)
    df
}

.checkSex <- function(sex) {
    if(! (is.character(sex) && !any(is.na(sex)) && all(sex %in% c("M", "F"))))
        stop("'sex' seems wrong (needs to be a character, without missing values, of 'M' and 'F'")
    sex
}

