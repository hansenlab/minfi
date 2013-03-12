library(minfiData)
library(digest)
Mset.raw <- preprocessRaw(RGsetEx)
Mset.swan <- preprocessSWAN(RGsetEx)
Mset.illumina <- preprocessSWAN(RGsetEx)

testDigests <- list(
    preprocessRaw = list(Meth = digest(getMeth(Mset.raw)),
    Unmeth = digest(getUnmeth(Mset.raw))),
    preprocessIllumina = list(Meth = digest(getMeth(Mset.illumina)),
    Unmeth = digest(getUnmeth(Mset.illumina))),
    preprocessSwan = list(Meth = digest(getMeth(Mset.swan)),
    Unmeth = digest(getUnmeth(Mset.swan)))
    )

save(testDigests, file = "../unitTests/testDigests.rda")
                         
