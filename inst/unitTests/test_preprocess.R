## test_subsetRGChannelSet <- function() {
##     stopifnot(require(minfiData))
##     exclude <- c("cg02004872", #  I-grn
##                  "cg00050873", # I-red
##                  "cg00035864") # II
##     rgSetSub <- subsetByLoci(RGsetEx, excludeLoci = exclude)
##     checkEquals(dim(preprocessRaw(rgSetSub)), dim(MsetEx[! rownames(MsetEx) %in% exclude,]))
## }
