test_sex <- function() {
    stopifnot(require(minfiData))
    data(MsetEx)
    gmSet <- mapToGenome(MsetEx)
    gmSet <- addSex(gmSet)
    checkEquals(gmSet$sex, gmSet$predictedSex)
}

