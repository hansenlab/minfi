library(minfiData)
library(digest)

GRset.funnorm <- preprocessFunnorm(RGsetEx)
digest_preprocessFunnorm <- list(M = minfi:::.digestMatrix(getM(GRset.funnorm)),
                                 CN = minfi:::.digestMatrix(getCN(GRset.funnorm)))
save(digest_preprocessFunnorm, file = "../unitTests/digest_preprocessFunnorm.rda")

sessionInfo()
rm(list = ls())

# ------------------------------------------------------------------------------
# DelayedArray tests
#

library(DelayedArray)

assays(RGsetEx) <- endoapply(assays(RGsetEx), DelayedArray)

GRset.funnorm <- preprocessFunnorm(RGsetEx)
digest_preprocessFunnorm <- list(M = minfi:::.digestMatrix(getM(GRset.funnorm)),
                                 CN = minfi:::.digestMatrix(getCN(GRset.funnorm)))
save(digest_preprocessFunnorm, file = "../unitTests/digest_preprocessFunnorm.DelayedArray.rda")

sessionInfo()
rm(list = ls())
