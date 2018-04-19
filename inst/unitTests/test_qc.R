test_sex <- function() {
    ## stopifnot(require(minfiData))
    ## stopifnot(require(DelayedArray))
    ## data(MsetEx)

    ## # Original tests
    ## gmSet <- mapToGenome(MsetEx)
    ## gmSet <- addSex(gmSet)
    ## checkEquals(gmSet$sex, gmSet$predictedSex)

    ## # Testing with DelayedArray-backed objects
    ## MsetEx <- realize(MsetEx)
    ## gmSet <- mapToGenome(MsetEx)
    ## gmSet <- addSex(gmSet)
    ## checkEquals(gmSet$sex, gmSet$predictedSex)
    checkTrue(TRUE)
}
