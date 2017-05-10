test_tcgaPipeline <- function() {
    stopifnot(require(minfiData))
    grSetEx <- tcgaPipeline(updateObject(RGsetEx))
    checkIdentical(getSnpBeta(grSetEx), getSnpBeta(updateObject(RGsetEx)))
    checkIdentical(preprocessMethod(grSetEx), 
                   c(mu.norm="Noob, dyeCorr=TRUE, dyeMethod=single")) 
}
