library(minfiData)
library(digest)

GMsetEx <- mapToGenome(MsetEx)
gr.cor <- createCorMatrix(GMsetEx, res=500*1000)
set.seed(456)
gr.ab <- extractAB(gr.cor)

digest_compartments <- list(cor.matrix = minfi:::.digestMatrix(gr.cor$cor.matrix),
                           pc = minfi:::.digestVector(gr.ab$pc, digits = 2))
save(digest_compartments, file = "../unitTests/digest_compartments.rda")

sessionInfo()                         
rm(list = ls())
