library(minfiData)
library(digest)

Mset.raw <- preprocessRaw(RGsetEx)
digest_preprocessRaw <- list(Meth = minfi:::.digestMatrix(getMeth(Mset.raw)),
                              Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.raw)))
save(digest_preprocessRaw, file = "../unitTests/digest_preprocessRaw.rda")

sessionInfo()                         
rm(list = ls())
