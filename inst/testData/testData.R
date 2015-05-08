library(minfiData)
library(digest)
Mset.raw <- preprocessRaw(RGsetEx)
Mset.illumina <- preprocessIllumina(RGsetEx)
set.seed(456)
Mset.swan <- preprocessSWAN(RGsetEx)
set.seed(456)
GRset.quantile <- preprocessQuantile(MsetEx)
set.seed(456)
Mset.noob <- preprocessNoob(RGsetEx)
set.seed(456)
GRset.funnorm <- preprocessFunnorm(RGsetEx)

gr.cor <- minfi:::createCorMatrix(MsetEx)
gr.ab <- minfi:::extractAB(gr.cor)

## save(Mset.raw, file = "Mset.raw.rda")
## save(Mset.illumina, file = "Mset.illumina.rda")
## save(Mset.swan, file = "Mset.swan.rda")
## save(Mset.quantile, file = "Mset.quantile.rda")


testDigests <- list(
    preprocessRaw = list(Meth = minfi:::.digestMatrix(getMeth(Mset.raw)),
      Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.raw))),
    preprocessIllumina = list(Meth = minfi:::.digestMatrix(getMeth(Mset.illumina)),
      Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.illumina))),
    preprocessSWAN = list(Meth = minfi:::.digestMatrix(getMeth(Mset.swan)),
      Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.swan))),
    preprocessQuantile = list(M = minfi:::.digestMatrix(getM(GRset.quantile)),
      CN = minfi:::.digestMatrix(getCN(GRset.quantile))),
    preprocessNoob = list(Meth = minfi:::.digestMatrix(getMeth(Mset.noob)),
      Unmeth = minfi:::.digestMatrix(getUnmeth(Mset.noob))),
    preprocessFunnorm = list(M = minfi:::.digestMatrix(getM(GRset.funnorm)),
      CN = minfi:::.digestMatrix(getCN(GRset.funnorm))),
    createCorMatrix = list(cor.matrix = minfi:::.digestMatrix(gr.cor$cor.matrix)),
    extractAB = list(pc = minfi:::.digestVector(gr.ab$pc))
)

save(testDigests, file = "../unitTests/testDigests.rda")

gc()
sessionInfo()                         

rm(list = ls())
