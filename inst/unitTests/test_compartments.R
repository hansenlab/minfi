test_compartments <- function() {
    ## stopifnot(require(minfiData))
    ## stopifnot(require(digest))
    ## stopifnot(require(DelayedArray))
    ## load(file.path(path.package("minfi"),
    ##                "unitTests",
    ##                "digest_compartments.rda"))

    ## # Original test
    ## GMsetEx <- mapToGenome(MsetEx)
    ## gr.cor <- createCorMatrix(GMsetEx, res = 500*1000)
    ## checkEquals(digest_compartments$cor.matrix,
    ##             minfi:::.digestMatrix(gr.cor$cor.matrix))
    ## set.seed(456)
    ## gr.ab <- extractAB(gr.cor)
    ## checkEquals(digest_compartments$pc,
    ##             minfi:::.digestVector(gr.ab$pc, digits = 2))

    ## # Testing with DelayedArray-backed objects
    ## MsetEx <- realize(MsetEx)
    ## GMsetEx <- mapToGenome(MsetEx)
    ## checkException(gr.cor <- createCorMatrix(GMsetEx, res = 500*1000),
    ##                silent = TRUE)
    ## # TODO: Uncomment once createCorMatrix() supports DelayedArray-backed
    ## #       minfi objects
    ## # gr.cor <- createCorMatrix(GMsetEx, res = 500*1000)
    ## # checkEquals(digest_compartments$cor.matrix,
    ## #             minfi:::.digestMatrix(gr.cor$cor.matrix))
    ## # set.seed(456)
    ## # gr.ab <- extractAB(gr.cor)
    ## # checkEquals(digest_compartments$pc,
    ## #             minfi:::.digestVector(gr.ab$pc, digits = 2))
    checkTrue(TRUE)
}
