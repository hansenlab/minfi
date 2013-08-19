library(minfiData)
library(digest)
Mset.raw <- preprocessRaw(RGsetEx)
Mset.illumina <- preprocessIllumina(RGsetEx)
set.seed(456)
Mset.swan <- preprocessSWAN(RGsetEx)
Mset.quantile <- preprocessQuantile(MsetEx)

## save(Mset.raw, file = "Mset.raw.rda")
## save(Mset.illumina, file = "Mset.illumina.rda")
## save(Mset.swan, file = "Mset.swan.rda")
save(Mset.quantile, file = "Mset.quantile.rda")


testDigests <- list(
    preprocessRaw = list(Meth = minfi:::.digestMatrix(getMeth(Mset.raw)),
      Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.raw))),
    preprocessIllumina = list(Meth = minfi:::.digestMatrix(getMeth(Mset.illumina)),
      Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.illumina))),
    preprocessSWAN = list(Meth = minfi:::.digestMatrix(getMeth(Mset.swan)),
      Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.swan))),
    preprocessQuantile = list(M = minfi:::.digestMatrix(getM(Mset.quantile)),
      CN = minfi:::.digestMatrix(getCN(Mset.quantile)))
    )

save(testDigests, file = "../unitTests/testDigests.rda")

gc()
sessionInfo()                         

rm(list = ls())
