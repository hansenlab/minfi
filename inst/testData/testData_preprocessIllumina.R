library(minfiData)
library(digest)

Mset.illumina <- preprocessIllumina(RGsetEx)
digest_preprocessIllumina <- list(Meth = minfi:::.digestMatrix(getMeth(Mset.illumina)),
                                      Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.illumina)))
save(digest_preprocessIllumina, file = "../unitTests/digest_preprocessIllumina.rda")

sessionInfo()
rm(list = ls())

# ------------------------------------------------------------------------------
# DelayedArray tests
#

library(DelayedArray)

assays(RGsetEx) <- endoapply(assays(RGsetEx), DelayedArray)

Mset.illumina <- preprocessIllumina(RGsetEx)
digest_preprocessIllumina <- list(Meth = minfi:::.digestMatrix(getMeth(Mset.illumina)),
                                  Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.illumina)))
save(digest_preprocessIllumina, file = "../unitTests/digest_preprocessIllumina.DelayedArray.rda")

sessionInfo()
rm(list = ls())
