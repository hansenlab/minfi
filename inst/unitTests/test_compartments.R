test_compartments <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_compartments.rda"))
    GMsetEx <- mapToGenome(MsetEx)
    gr.cor <- createCorMatrix(GMsetEx, res=500*1000)
    checkEquals(digest_compartments$cor.matrix,
                minfi:::.digestMatrix(gr.cor$cor.matrix))
    set.seed(456)
    gr.ab <- extractAB(gr.cor)
    checkEquals(digest_compartments$pc,
                minfi:::.digestVector(gr.ab$pc, digits = 2))
}

# ------------------------------------------------------------------------------
# DelayedArray tests
#

library(DelayedArray)

test_compartments_with_DelayedArray <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("minfi"), "unitTests", "digest_compartments.DelayedArray.rda"))
    assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)
    GMsetEx <- mapToGenome(MsetEx)
    gr.cor <- createCorMatrix(GMsetEx, res=500*1000)
    checkEquals(digest_compartments$cor.matrix,
                minfi:::.digestMatrix(gr.cor$cor.matrix))
    set.seed(456)
    gr.ab <- extractAB(gr.cor)
    checkEquals(digest_compartments$pc,
                minfi:::.digestVector(gr.ab$pc, digits = 2))
}
