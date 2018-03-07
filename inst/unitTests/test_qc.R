test_sex <- function() {
    stopifnot(require(minfiData))
    data(MsetEx)
    gmSet <- mapToGenome(MsetEx)
    gmSet <- addSex(gmSet)
    checkEquals(gmSet$sex, gmSet$predictedSex)
}

# ------------------------------------------------------------------------------
# DelayedArray tests
#

library(DelayedArray)

test_sex_with_DelayedArray <- function() {
    stopifnot(require(minfiData))
    data(MsetEx)

    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)

    gmSet <- mapToGenome(MsetEx)
    gmSet <- addSex(gmSet)
    checkEquals(gmSet$sex, gmSet$predictedSex)
}
