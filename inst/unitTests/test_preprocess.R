test_preprocessRaw <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "testDigests.rda"))
    Mset <- preprocessRaw(RGsetEx)
    checkEquals(testDigests$preprocessRaw$Meth,
                digest(getMeth(Mset)))
    checkEquals(testDigests$preprocessRaw$Unmeth,
                digest(getUnmeth(Mset)))
}

test_preprocessIllumina <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "testDigests.rda"))
    Mset <- preprocessIllumina(RGsetEx)
    checkEquals(testDigests$preprocessIllumina$Meth,
                digest(getMeth(Mset)))
    checkEquals(testDigests$preprocessIllumina$Unmeth,
                digest(getUnmeth(Mset)))
}

test_preprocessSWAN <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "testDigests.rda"))
    set.seed(456)
    Mset <- preprocessSWAN(RGsetEx)
    checkEquals(testDigests$preprocessSWAN$Meth,
                digest(getMeth(Mset)))
    checkEquals(testDigests$preprocessSWAN$Unmeth,
                digest(getUnmeth(Mset)))
}



