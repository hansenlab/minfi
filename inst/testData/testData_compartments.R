library(minfiData)
library(digest)

gr.cor <- minfi:::createCorMatrix(MsetEx, res=500*1000)
gr.ab <- minfi:::extractAB(gr.cor)

digest_compartments <- list(cor.matrix = minfi:::.digestMatrix(gr.cor$cor.matrix),
                           pc = minfi:::.digestVector(gr.ab$pc, digits = 3))
save(digest_compartments, file = "../unitTests/digest_compartments.rda")

sessionInfo()                         
rm(list = ls())
