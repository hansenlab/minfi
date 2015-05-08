library(minfiData)
library(digest)

set.seed(456)
Mset.swan <- preprocessSWAN(RGsetEx)
digest_preprocessSWAN <- list(Meth = minfi:::.digestMatrix(getMeth(Mset.swan)),
                                  Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.swan)))
save(digest_preprocessSWAN, file = "../unitTests/digest_preprocessSWAN.rda")

sessionInfo()                         
rm(list = ls())
