library(minfiData)
library(digest)

GRset.quantile <- preprocessQuantile(MsetEx)
digest_preprocessQuantile <- list(M = minfi:::.digestMatrix(getM(GRset.quantile)),
                                  CN = minfi:::.digestMatrix(getCN(GRset.quantile)))
save(digest_preprocessQuantile, file = "../unitTests/digest_preprocessQuantile.rda")

sessionInfo()
rm(list = ls())

# ------------------------------------------------------------------------------
# DelayedArray tests
#

library(DelayedArray)

assays(MsetEx) <- endoapply(assays(MsetEx), DelayedArray)

GRset.quantile <- preprocessQuantile(MsetEx)
digest_preprocessQuantile <- list(M = minfi:::.digestMatrix(getM(GRset.quantile)),
                                  CN = minfi:::.digestMatrix(getCN(GRset.quantile)))
save(digest_preprocessQuantile, file = "../unitTests/digest_preprocessQuantile.DelayedArray.rda")

sessionInfo()
rm(list = ls())
