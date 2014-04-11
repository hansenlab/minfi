test_preprocessRaw <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "testDigests.rda"))
    Mset <- preprocessRaw(RGsetEx)
    checkEquals(testDigests$preprocessRaw$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(testDigests$preprocessRaw$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))
}

test_preprocessIllumina <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "testDigests.rda"))
    Mset <- preprocessIllumina(RGsetEx)
    checkEquals(testDigests$preprocessIllumina$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(testDigests$preprocessIllumina$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))
}

test_preprocessSWAN <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "testDigests.rda"))
    set.seed(456)
    Mset <- preprocessSWAN(RGsetEx)
    checkEquals(testDigests$preprocessSWAN$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(testDigests$preprocessSWAN$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))
}

test_preprocessQuantile <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "testDigests.rda"))
    GMset <- preprocessQuantile(MsetEx)
    checkEquals(testDigests$preprocessQuantile$M,
                minfi:::.digestMatrix(getM(GMset)))
    checkEquals(testDigests$preprocessQuantile$CN,
                minfi:::.digestMatrix(getCN(GMset)))
}



