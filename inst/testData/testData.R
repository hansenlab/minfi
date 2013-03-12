library(minfiData)
library(digest)
Mset.raw <- preprocessRaw(RGsetEx)
Mset.swan <- preprocessSWAN(RGsetEx)
set.seed(456)
Mset.illumina <- preprocessIllumina(RGsetEx)

testDigests <- list(
    preprocessRaw = list(Meth = digest(getMeth(Mset.raw)),
    Unmeth = digest(getUnmeth(Mset.raw))),
    preprocessIllumina = list(Meth = digest(getMeth(Mset.illumina)),
    Unmeth = digest(getUnmeth(Mset.illumina))),
    preprocessSWAN = list(Meth = digest(getMeth(Mset.swan)),
    Unmeth = digest(getUnmeth(Mset.swan)))
    )

save(testDigests, file = "../unitTests/testDigests.rda")

sessionInfo()                         
