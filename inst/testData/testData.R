library(minfiData)
library(digest)
Mset.raw <- preprocessRaw(RGsetEx)
Mset.illumina <- preprocessIllumina(RGsetEx)
set.seed(456)
Mset.swan <- preprocessSWAN(RGsetEx)

testDigests <- list(
    preprocessRaw = list(Meth = minfi:::.digestMatrix(getMeth(Mset.raw)),
    Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.raw))),
    preprocessIllumina = list(Meth = minfi:::.digestMatrix(getMeth(Mset.illumina)),
    Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.illumina))),
    preprocessSWAN = list(Meth = minfi:::.digestMatrix(getMeth(Mset.swan)),
    Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.swan)))
    )

save(testDigests, file = "../unitTests/testDigests.rda")

gc()
sessionInfo()                         
