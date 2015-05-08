library(minfiData)
library(digest)

GRset.quantile <- preprocessQuantile(MsetEx)
digest_preprocessQuantile <- list(M = minfi:::.digestMatrix(getM(GRset.quantile)),
                                  CN = minfi:::.digestMatrix(getCN(GRset.quantile)))
save(digest_preprocessQuantile, file = "../unitTests/digest_preprocessQuantile.rda")

sessionInfo()                         
rm(list = ls())
