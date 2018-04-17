# Internal functions -----------------------------------------------------------

.getSex <- function(CN = NULL, xIndex = NULL, yIndex = NULL, cutoff = -2) {
    if (is.null(CN) || is.null(xIndex) || is.null(yIndex)) {
        stop("must provide CN, xIndex, and yIndex")
    }
    # TODO: This does not handle only females or only males. This ought to be
    #       handled by the 'centers' (see below) being too close together
    xMed <- colMedians(CN, rows = xIndex, na.rm = TRUE)
    yMed <- colMedians(CN, rows = yIndex, na.rm = TRUE)
    dd <- yMed - xMed
    k <- kmeans(dd, centers = c(min(dd), max(dd)))

    sex0 <- ifelse(dd < cutoff, "F", "M")
    sex0 <- .checkSex(sex0)
    sex1 <- ifelse(k$cluster == which.min(k$centers), "F", "M")
    sex1 <- .checkSex(sex1)

    if (!identical(sex0, sex1)) {
        warning("An inconsistency was encountered while determining sex. One ",
                "possibility is that only one sex is present. We recommend ",
                "further checks, for example with the plotSex function.")
    }
    df <- DataFrame(xMed = xMed, yMed = yMed, predictedSex = sex0)
    rownames(df) <- colnames(CN)
    df
}

.checkSex <- function(sex) {
    if (!(is.character(sex) && !any(is.na(sex)) &&
          all(sex %in% c("M", "F")))) {
        stop("'sex' seems wrong (needs to be a character, without ",
             "missing values, of 'M' and 'F'")
    }
    sex
}

# Exported functions -----------------------------------------------------------

getSex <- function(object = NULL, cutoff = -2){
    .isGenomicOrStop(object)
    if (is(object, "GenomicMethylSet")) CN <- getCN(object)
    if (is(object, "GenomicRatioSet")) CN <- getCN(object)
    # TODO: Add test for logarithmic scale or non-log scale
    xIndex <- which(seqnames(object) == "chrX")
    yIndex <- which(seqnames(object) == "chrY")
    .getSex(
        CN = CN,
        xIndex = xIndex,
        yIndex = yIndex,
        cutoff = cutoff)
}

addSex <- function(object, sex = NULL) {
    if (is.null(sex)) {
        sex <- getSex(object)
    }
    if (is(sex, "DataFrame")) {
        stopifnot(all(c("predictedSex", "xMed", "yMed") %in% colnames(sex)))
    }
    sex$predictedSex <- .checkSex(sex$predictedSex)
    .pDataAdd(object, sex)
}

plotSex <- function(object, id = NULL) {
    stopifnot(all(c("predictedSex", "xMed", "yMed") %in%
                      colnames(colData(object))))
    if (is.null(id)) id <- seq_along(object$predictedSex)
    if (length(id) != length(object$predictedSex)) {
        stop("id length must match number of samples.")
    }
    plot(
        x = object$xMed,
        y = object$yMed,
        type = "n",
        xlab = "X chr, median total intensity (log2)",
        ylab = "Y chr, median total intensity (log2)")
    text(
        x = object$xMed,
        y = object$yMed,
        labels = id,
        col = ifelse(object$predictedSex == "M", "deepskyblue", "deeppink3"))
    legend(
        "bottomleft",
        c("M", "F"),
        col = c("deepskyblue", "deeppink3"),
        pch = 16)
}
