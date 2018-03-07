library(minfiData)
library(digest)

Mset.noob <- preprocessNoob(RGsetEx)
digest_preprocessNoob <- list(Meth = minfi:::.digestMatrix(getMeth(Mset.noob)),
                              Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.noob)))
save(digest_preprocessNoob, file = "../unitTests/digest_preprocessNoob.rda")

sessionInfo()
rm(list = ls())

# ------------------------------------------------------------------------------
# DelayedArray tests
#

library(DelayedArray)

assays(RGsetEx) <- endoapply(assays(RGsetEx), DelayedArray)
Mset.noob <- preprocessNoob(RGsetEx)
digest_preprocessNoob <- list(Meth = minfi:::.digestMatrix(getMeth(Mset.noob)),
                              Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.noob)))
save(digest_preprocessNoob, file = "../unitTests/digest_preprocessNoob.DelayedArray.rda")

sessionInfo()
rm(list = ls())
