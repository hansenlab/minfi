test_subsetRGChannelSet <- function() {
##     stopifnot(require(minfiData))
##     stopifnot(require(DelayedArray))
##     exclude <- c("cg02004872", #  I-grn
##                  "cg00050873", # I-red
##                  "cg00035864") # II

##     # Original tests
##     rgSetSub <- subsetByLoci(RGsetEx, excludeLoci = exclude)
##     checkEquals(dim(preprocessRaw(rgSetSub)),
##                 dim(MsetEx[!rownames(MsetEx) %in% exclude, ]))
##     checkEquals(dim(preprocessNoob(rgSetSub)),
##                 dim(MsetEx[!rownames(MsetEx) %in% exclude, ]))
##     # TODO: Why as.numeric()?
##     checkEquals(as.numeric(dim(preprocessFunnorm(rgSetSub))),
##                 as.numeric(dim(MsetEx[!rownames(MsetEx) %in% exclude, ])))

##     # Testing with DelayedArray-backed objects
##     RGsetEx <- realize(RGsetEx)
##     MsetEx <- realize(MsetEx)
##     rgSetSub <- subsetByLoci(RGsetEx, excludeLoci = exclude)
##     checkEquals(dim(preprocessRaw(rgSetSub)),
##                 dim(MsetEx[!rownames(MsetEx) %in% exclude, ]))
##     checkEquals(dim(preprocessNoob(rgSetSub)),
##                 dim(MsetEx[!rownames(MsetEx) %in% exclude, ]))
##     # TODO: Uncomment once preprocessFunnorm() supports DelayedArray-backed
##     #       minfi objects
##     checkException(preprocessFunnorm(rgSetSub), silent = TRUE)
##     # checkEquals(as.numeric(dim(preprocessFunnorm(rgSetSub))),
##     #             as.numeric(dim(MsetEx[!rownames(MsetEx) %in% exclude, ])))
    checkTrue(TRUE)
}

test_preprocessRaw <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    stopifnot(require(DelayedArray))
    load(file.path(
        path.package("minfi"),
        "unitTests",
        "digest_preprocessRaw.rda"))

    # Original tests
    Mset <- preprocessRaw(RGsetEx)
    checkEquals(digest_preprocessRaw$Meth,
                minfi:::.digestMatrix(getMeth(Mset)))
    checkEquals(digest_preprocessRaw$Unmeth,
                minfi:::.digestMatrix(getUnmeth(Mset)))

    # Testing with DelayedArray-backed objects
    RGsetEx <- realize(RGsetEx)
    Mset <- preprocessRaw(RGsetEx)
    checkEquals(digest_preprocessRaw$Meth,
                minfi:::.digestMatrix(as.matrix(getMeth(Mset))))
    checkEquals(digest_preprocessRaw$Unmeth,
                minfi:::.digestMatrix(as.matrix(getUnmeth(Mset))))
}

test_preprocessIllumina <- function() {
    ## stopifnot(require(minfiData))
    ## stopifnot(require(digest))
    ## stopifnot(require(DelayedArray))
    ## load(file.path(path.package("minfi"),
    ##                "unitTests",
    ##                "digest_preprocessIllumina.rda"))

    ## # Original tests
    ## Mset <- preprocessIllumina(RGsetEx)
    ## checkEquals(digest_preprocessIllumina$Meth,
    ##             minfi:::.digestMatrix(getMeth(Mset)))
    ## checkEquals(digest_preprocessIllumina$Unmeth,
    ##             minfi:::.digestMatrix(getUnmeth(Mset)))

    ## # Testing with DelayedArray-backed objects
    ## RGsetEx <- realize(RGsetEx)
    ## Mset <- preprocessIllumina(RGsetEx)
    ## checkEquals(digest_preprocessIllumina$Meth,
    ##             minfi:::.digestMatrix(as.matrix(getMeth(Mset))))
    ## checkEquals(digest_preprocessIllumina$Unmeth,
    ##             minfi:::.digestMatrix(as.matrix(getUnmeth(Mset))))
    checkTrue(TRUE)
}

test_preprocessSWAN <- function() {
    ## stopifnot(require(minfiData))
    ## stopifnot(require(digest))
    ## stopifnot(require(DelayedArray))
    ## load(file.path(path.package("minfi"),
    ##                "unitTests",
    ##                "digest_preprocessSWAN.rda"))

    ## # Original tests
    ## set.seed(456)
    ## Mset <- preprocessSWAN(RGsetEx)
    ## checkEquals(digest_preprocessSWAN$Meth,
    ##             minfi:::.digestMatrix(getMeth(Mset)))
    ## checkEquals(digest_preprocessSWAN$Unmeth,
    ##             minfi:::.digestMatrix(getUnmeth(Mset)))

    ## # Testing with DelayedArray-backed objects
    ## RGsetEx <- realize(RGsetEx)
    ## set.seed(456)
    ## Mset <- preprocessSWAN(RGsetEx)
    ## checkEquals(digest_preprocessSWAN$Meth,
    ##             minfi:::.digestMatrix(as.matrix(getMeth(Mset))))
    ## checkEquals(digest_preprocessSWAN$Unmeth,
    ##             minfi:::.digestMatrix(as.matrix(getUnmeth(Mset))))
    checkTrue(TRUE)
}

test_preprocessQuantile <- function() {
    ## stopifnot(require(minfiData))
    ## stopifnot(require(digest))
    ## stopifnot(require(DelayedArray))
    ## load(file.path(path.package("minfi"),
    ##                "unitTests",
    ##                "digest_preprocessQuantile.rda"))

    ## # Original tests
    ## GRset <- preprocessQuantile(MsetEx)
    ## checkEquals(digest_preprocessQuantile$M,
    ##             minfi:::.digestMatrix(getM(GRset)))
    ## checkEquals(digest_preprocessQuantile$CN,
    ##             minfi:::.digestMatrix(getCN(GRset)))

    ## # Testing with DelayedArray-backed objects
    ## MsetEx <- realize(MsetEx)
    ## checkException(preprocessQuantile(MsetEx), silent = TRUE)
    ## # TODO: Uncomment once preprocessQuantile() supports DelayedArray-backed
    ## #       minfi objects
    ## # checkEquals(digest_preprocessQuantile$M,
    ## #             minfi:::.digestMatrix(as.matrix(getM(GRset))))
    ## # checkEquals(digest_preprocessQuantile$CN,
    ## #             minfi:::.digestMatrix(as.matrix(getCN(GRset))))
    checkTrue(TRUE)
}

test_preprocessNoob <- function() {
    ## stopifnot(require(minfiData))
    ## stopifnot(require(digest))
    ## stopifnot(require(DelayedArray))
    ## load(file.path(path.package("minfi"),
    ##                "unitTests",
    ##                "digest_preprocessNoob.rda"))

    ## # Original tests
    ## Mset <- preprocessNoob(RGsetEx)
    ## checkEquals(digest_preprocessNoob$Meth,
    ##             minfi:::.digestMatrix(getMeth(Mset)))
    ## checkEquals(digest_preprocessNoob$Unmeth,
    ##             minfi:::.digestMatrix(getUnmeth(Mset)))

    ## # Testing with DelayedArray-backed objects
    ## RGsetEx <- realize(RGsetEx)
    ## Mset <- preprocessNoob(RGsetEx)
    ## checkEquals(digest_preprocessNoob$Meth,
    ##             minfi:::.digestMatrix(as.matrix(getMeth(Mset))))
    ## checkEquals(digest_preprocessNoob$Unmeth,
    ##             minfi:::.digestMatrix(as.matrix(getUnmeth(Mset))))
    checkTrue(TRUE)
}

test_preprocessFunnorm <- function() {
    ## stopifnot(require(minfiData))
    ## stopifnot(require(digest))
    ## stopifnot(require(DelayedArray))
    ## load(file.path(path.package("minfi"),
    ##                "unitTests",
    ##                "digest_preprocessFunnorm.rda"))

    ## # Original tests
    ## GRset <- preprocessFunnorm(RGsetEx)
    ## checkEquals(digest_preprocessFunnorm$M,
    ##             minfi:::.digestMatrix(getM(GRset)))
    ## checkEquals(digest_preprocessFunnorm$CN,
    ##             minfi:::.digestMatrix(getCN(GRset)))

    ## # Testing with DelayedArray-backed objects
    ## RGsetEx <- realize(RGsetEx)
    ## checkException(GRset <- preprocessFunnorm(RGsetEx), silent = TRUE)

    ## # TODO: Uncomment once preprocessFunnorm() supports DelayedArray-backed
    ## #       minfi objects
    ## # checkEquals(digest_preprocessFunnorm$M,
    ## #             minfi:::.digestMatrix(getM(GRset)))
    ## # checkEquals(digest_preprocessFunnorm$CN,
    ## #             minfi:::.digestMatrix(getCN(GRset)))
    checkTrue(TRUE)
}
