library(minfiData)
library(digest)

gr.cor <- createCorMatrix(MsetEx, res=500*1000)
gr.ab <- extractAB(gr.cor)

digest_compartents <- list(cor.matrix = minfi:::.digestMatrix(gr.cor$cor.matrix),
                           pc = minfi:::.digestVec(gr.ab$pc))
save(digest_compartments, file = "../unitTests/digest_compartments.rda")

sessionInfo()                         
rm(list = ls())
