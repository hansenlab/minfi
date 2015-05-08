test_preprocessRaw <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessRaw.rda"))
    Mset <- preprocessRaw(RGsetEx)
    checkEquals(digest_preprocessRaw$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(digest_preprocessRaw$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))
}

test_preprocessIllumina <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessIllumina.rda"))
    Mset <- preprocessIllumina(RGsetEx)
    checkEquals(digest_preprocessIllumina$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(digest_preprocessIllumina$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))
}

test_preprocessSWAN <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessSWAN.rda"))
    set.seed(456)
    Mset <- preprocessSWAN(RGsetEx)
    checkEquals(digest_preprocessSWAN$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(digest_preprocessSWAN$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))
}

test_preprocessQuantile <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessQuantile.rda"))
    GRset <- preprocessQuantile(MsetEx)
    checkEquals(digest_preprocessQuantile$M,
                minfi:::.digestMatrix(getM(GRset)))
    checkEquals(digest_preprocessQuantile$CN,
                minfi:::.digestMatrix(getCN(GRset)))
}

test_preprocessNoob <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessNoob.rda"))
    Mset <- preprocessNoob(RGsetEx)
    checkEquals(digest_preprocessNoob$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(digest_preprocessNoob$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))
}

test_preprocessFunnorm <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessFunnorm.rda"))
    GRset <- preprocessFunnorm(RGsetEx)
    checkEquals(digest_preprocessFunnorm$M,
                minfi:::.digestMatrix(getM(GRset)))
    checkEquals(digest_preprocessFunnorm$CN,
                minfi:::.digestMatrix(getCN(GRset)))
}

