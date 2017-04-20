test_noobPipeline <- function() {
    stopifnot(require(minfiData))
    grSetEx <- noobPipeline(updateObject(RGsetEx))
    checkIdentical(getSNPs(grSetEx), getSNPs(updateObject(RGsetEx)))
    checkIdentical(preprocessMethod(grSetEx), 
                   c(mu.norm="Noob, dyeCorr=TRUE, dyeMethod=single")) 
}
