test_subsetRGChannelSet <- function() {
    stopifnot(require(minfiData))
    exclude <- c("cg02004872", #  I-grn
                 "cg00050873", # I-red
                 "cg00035864") # II
    rgSetSub <- subsetByLoci(RGsetEx, excludeLoci = exclude)
    checkEquals(dim(preprocessRaw(rgSetSub)), dim(MsetEx[! rownames(MsetEx) %in% exclude,]))
    checkEquals(dim(preprocessNoob(rgSetSub)), dim(MsetEx[! rownames(MsetEx) %in% exclude,]))
    checkEquals(as.numeric(dim(preprocessFunnorm(rgSetSub))), as.numeric(dim(MsetEx[! rownames(MsetEx) %in% exclude,])))
}

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

# ------------------------------------------------------------------------------
# DelayedArray tests
#

library(DelayedArray)

test_subsetRGChannelSet_with_DelayedArray <- function() {
    stopifnot(require(minfiData))
    exclude<- c("cg02004872", #  I-grn
                "cg00050873", # I-red
                "cg00035864") # II

    assays(RGsetEx) <- endoapply(assays(RGsetEx, DelayedArray))

    rgSetSub <- subsetByLoci(RGsetEx, excludeLoci = exclude)
    checkEquals(dim(preprocessRaw(rgSetSub)), dim(MsetEx[! rownames(MsetEx) %in% exclude,]))
    checkEquals(dim(preprocessNoob(rgSetSub)), dim(MsetEx[! rownames(MsetEx) %in% exclude,]))
    checkEquals(as.numeric(dim(preprocessFunnorm(rgSetSub))), as.numeric(dim(MsetEx[! rownames(MsetEx) %in% exclude,])))
}

test_preprocessRaw_with_DelayedArray <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessRaw.DelayedArray.rda"))

    assays(RGsetEx) <- endoapply(assays(RGsetEx, DelayedArray))

    Mset <- preprocessRaw(RGsetEx)
    checkEquals(digest_preprocessRaw$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(digest_preprocessRaw$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))
}

test_preprocessIllumina_with_DelayedArray <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessIllumina.DelayedArray.rda"))

    assays(RGsetEx) <- endoapply(assays(RGsetEx, DelayedArray))

    Mset <- preprocessIllumina(RGsetEx)
    checkEquals(digest_preprocessIllumina$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(digest_preprocessIllumina$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))
}

test_preprocessSWAN_with_DelayedArray <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessSWAN.DelayedArray.rda"))
    set.seed(456)

    assays(RGsetEx) <- endoapply(assays(RGsetEx, DelayedArray))

    Mset <- preprocessSWAN(RGsetEx)
    checkEquals(digest_preprocessSWAN$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(digest_preprocessSWAN$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))
}

test_preprocessQuantile_with_DelayedArray <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessQuantile.DelayedArray.rda"))

    assays(MsetEx) <- endoapply(assays(MsetEx, DelayedArray))

    GRset <- preprocessQuantile(MsetEx)
    checkEquals(digest_preprocessQuantile$M,
                minfi:::.digestMatrix(getM(GRset)))
    checkEquals(digest_preprocessQuantile$CN,
                minfi:::.digestMatrix(getCN(GRset)))
}

test_preprocessNoob_with_DelayedArray <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessNoob.DelayedArray.rda"))

    assays(RGsetEx) <- endoapply(assays(RGsetEx, DelayedArray))

    Mset <- preprocessNoob(RGsetEx)
    checkEquals(digest_preprocessNoob$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(digest_preprocessNoob$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))
}

test_preprocessFunnorm_with_DelayedArray <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_preprocessFunnorm.DelayedArray.rda"))

    assays(RGsetEx) <- endoapply(assays(RGsetEx, DelayedArray))

    GRset <- preprocessFunnorm(RGsetEx)
    checkEquals(digest_preprocessFunnorm$M,
                minfi:::.digestMatrix(getM(GRset)))
    checkEquals(digest_preprocessFunnorm$CN,
                minfi:::.digestMatrix(getCN(GRset)))
}
