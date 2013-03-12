test_preprocessRaw <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load("testDigests.rda")
    Mset <- preprocessRaw(RGsetEx)
    checkEquals(testDigests$preprocessRaw$Meth,
                digest(getMeth(Mset)))
    checkEquals(testDigests$preprocessRaw$Unmeth,
                digest(getUnmeth(Mset)))
}



