## This script compares the output of Genome Studio to preprocess.Illumina

require(minfi)
require(matrixStats)

basepath <- "/thumper2/feinbergLab/core/arrays/illumina/IL002"
sheet <- read.450k.sheet(basepath, pattern = "IL002_v2.csv$")
RGset <- read.450k.exp(targets = sheet, verbose = TRUE)
ISet.raw <- preprocessIllumina(RGset, bg.correct = FALSE, normalize = "no")
ISet.bg <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "no")
ISet.norm <- preprocessIllumina(RGset, bg.correct = FALSE, normalize = "controls",
                                reference = 8)
ISet.norm.bg <- preprocessIllumina(RGset, bg.correct = TRUE,
                                   normalize = "controls", reference = 8)

od <- c(1,7,8,2,14,9, 3,21,10,15,16,22,
        4,11,23,17,12,24,18,19,5,6,13,20)
basepath <- "/thumper2/feinbergLab/personal/khansen/450k/data"
GS.raw <- minfi:::read.GenomeStudio(file.path(basepath,
                                      "FinalReport_noBack_noNorm.txt"))
odN <- colnames(GS.raw$beta)[od]
GS.raw <- lapply(GS.raw, function(xx) as.matrix(xx[featureNames(ISet.raw), odN]))

GS.bg <- minfi:::read.GenomeStudio(file.path(basepath,
                                     "FinalReport_yesBack_noNorm.txt"))
GS.bg <- lapply(GS.bg, function(xx) as.matrix(xx[featureNames(ISet.bg), odN]))

GS.norm <- minfi:::read.GenomeStudio(file.path(basepath,
                                       "FinalReport_noBack_yesNorm.txt"))
GS.norm <- lapply(GS.norm, function(xx) as.matrix(xx[featureNames(ISet.norm), odN]))

GS.norm.bg <- minfi:::read.GenomeStudio(file.path(basepath,
                                          "FinalReport_yesBack_yesNorm.txt"))
GS.norm.bg <- lapply(GS.norm.bg, function(xx) as.matrix(xx[featureNames(ISet.norm.bg), odN]))

## Now for comparisons

cmax <- colMaxs(abs(getMeth(ISet.raw) - GS.raw$SignalB))
round(cmax, 1)
cmax <- colMaxs(abs(getUnmeth(ISet.raw) - GS.raw$SignalA))
round(cmax, 1)
cmax <- colMaxs(abs(getBeta(ISet.raw, type = "Illumina") - GS.raw$beta), na.rm = TRUE)
round(cmax, 3)

cmax <- colMaxs(abs(getMeth(ISet.bg) - GS.bg$SignalB))
round(cmax, 1)
cmax <- colMaxs(abs(getUnmeth(ISet.bg) - GS.bg$SignalA))
round(cmax, 1)
cmax <- colMaxs(abs(getBeta(ISet.bg, type = "Illumina") - GS.bg$beta), na.rm = TRUE)
round(cmax, 3)

cmax <- colMaxs(abs(getMeth(ISet.norm) - GS.norm$SignalB))
round(cmax, 2)
cmax <- colMaxs(abs(getUnmeth(ISet.norm) - GS.norm$SignalA))
round(cmax, 2)
cmax <- colMaxs(abs(getBeta(ISet.norm, type = "Illumina") - GS.norm$beta), na.rm = TRUE)
round(cmax, 3)

cmax <- colMaxs(abs(getMeth(ISet.norm.bg) - GS.norm.bg$SignalB))
round(cmax, 2)
cmax <- colMaxs(abs(getUnmeth(ISet.norm.bg) - GS.norm.bg$SignalA))
round(cmax, 2)
cmax <- colMaxs(abs(getBeta(ISet.norm.bg, type = "Illumina") - GS.norm.bg$beta), na.rm = TRUE)
round(cmax, 3)

