# R-2.13
setwd("~maryee/projects/infinium450k")
library(minfi)
#library(SQN)
x<-ls();for (file in list.files("functions", full.names=TRUE)) source(file);setdiff(ls(),x)
data(IlluminaHumanMethylation450kmanifest)
#library(Biostrings)

#manifestFile <- "/thumper2/feinbergLab/personal/khansen/450k/data/HumanMethylation450_15017482_v.1.1_ForExcel.csv"


################################################
## Effect of normalization options: MDS plots ##
################################################
## IL003
targets <- read.csv("data/IL003_detailed.csv", stringsAsFactors = FALSE, skip=7)
targets <- subset(targets, GrossAnatomicalSite %in% c("cerebellum", "cerebrum", "colon", "lymphocytes"))
load("rdas/rgset_il003.rda")

b <- getBeta(preprocess.raw(RGset))
bBg <- getBeta(preprocess.illumina(RGset, bg.correct=TRUE, normalize="no"))
bControl <- getBeta(preprocess.illumina(RGset, bg.correct=FALSE, normalize="controls", reference=1))
bBgControl <- getBeta(preprocess.illumina(RGset, bg.correct=TRUE, normalize="controls", reference=1))
bSqn <- getBeta(preprocess.raw(sqn(RGset)))
bBgSqn <- getBeta(preprocess.raw(sqn(bg.correct.illumina(RGset))))

## MDS plots
par(mfcol=c(2,3))
mdsPlot(b, n=1000, sampleGroups=targets$GrossAnatomicalSite, main="Raw\nBeta")
mdsPlot(bBg, n=1000, sampleGroups=targets$GrossAnatomicalSite, main="Background subtracted\nBeta")
mdsPlot(bControl, n=1000, sampleGroups=targets$GrossAnatomicalSite, main="Illumina control norm\nBeta")
mdsPlot(bBgControl, n=1000, sampleGroups=targets$GrossAnatomicalSite, main="Background + Illumina control norm\nBeta")
mdsPlot(bSqn, n=1000, sampleGroups=targets$GrossAnatomicalSite, main="SQN\nBeta")
mdsPlot(bBgSqn, n=1000, sampleGroups=targets$GrossAnatomicalSite, main="Background + SQN\nBeta")



##########################################
## Normalizing to a bad reference array ##
##########################################
## IL003
targets <- read.csv("data/IL003_detailed.csv", stringsAsFactors = FALSE, skip=7)
targets <- subset(targets, GrossAnatomicalSite %in% c("cerebellum", "cerebrum", "colon", "lymphocytes"))
load("rdas/rgset_il003.rda")

# Identify a bad sample using overall beta distribution
#mdsPlot(RGset, n=1000, sampleGroups=targets$GrossAnatomicalSite, sampleNames=targets$Sample_Name)
#badIdx <- which(targets$Sample_Name=="2993_autism_brain")

# Identify a bad sample using red channel controls
r <- getRed(RGset)
ctrlAddress <- getControlAddress(RGset, controlType=c("NORM_A", "NORM_T"))
y <- log2(r[ctrlAddress,])
#densityPlot(y, sampleGroups=targets$GrossAnatomicalSite)
badIdx1 <- which.min(colMeans(y))
badIdx2 <- which.max(colMeans(y))

b <- getBeta(preprocess.raw(RGset))
bControlBadRef1 <- getBeta(preprocess.illumina(RGset, bg.correct=FALSE, normalize="controls", reference=badIdx1))
bControlBadRef2 <- getBeta(preprocess.illumina(RGset, bg.correct=FALSE, normalize="controls", reference=badIdx2))

par(mfrow=c(2,2))
mdsPlot(b, n=1000, sampleGroups=targets$GrossAnatomicalSite)
mdsPlot(bControlBadRef1, n=1000, sampleGroups=targets$GrossAnatomicalSite)
mdsPlot(bControlBadRef2, n=1000, sampleGroups=targets$GrossAnatomicalSite)
