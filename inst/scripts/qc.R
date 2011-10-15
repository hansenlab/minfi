# R-2.13
setwd("~maryee/projects/infinium450k")
library(minfi)
for (file in list.files("functions", full.names=TRUE)) source(file)
data(IlluminaHumanMethylation450kmanifest)
#library(Biostrings)

###############
## Functions ##
###############
# None defined here


#############################
## replication_IL002_IL003 ##
#############################
targets <- read.csv("data/replication_IL002_IL003.csv", stringsAsFactors = FALSE, skip = 7)
rda <- "rdas/rgset_replication_il002_il003.rda"
load(rda)
sampleNames <- paste(targets$Experiment, targets$Sample_Name, targets$Sample_Well)

## Density plots
densityPlot(RGset, sampleGroups=targets$Sample_Group, main="Beta", xlab="Beta")

## Beta beanplots
par(oma=c(2,10,1,1))
densityBeanPlot(RGset, sampleGroups=targets$Sample_Group, sampleNames=sampleNames)

## Control probe stripplots
#par(ask=TRUE)
controlStripPlot(RGset, controls="BISULFITE CONVERSION II", sampleNames=sampleNames)

## Full qcReport
qcReport(RGset, sampleNames=sampleNames, sampleGroups=targets$Sample_Group, pdf="figs/qcReport_replication_IL002_IL003.pdf")

## MDS plots
mdsPlot(RGset, n=1000, sampleGroups=targets$Sample_Group)
mdsPlot(RGset, n=1000, sampleGroups=targets$Sample_Group, sampleNames=sampleNames)



###########
## IL003 ##
###########
targets <- read.csv("data/IL003_detailed.csv", stringsAsFactors = FALSE)
targets <- subset(targets, GrossAnatomicalSite %in% c("cerebellum", "cerebrum", "colon", "lymphocytes"))
rda <- "rdas/rgset_il003.rda"
load(rda)
sampleNames <- paste(targets$Sample_Name, targets$Sample_Well)
qcReport(RGset, sampleNames=sampleNames, sampleGroups=targets$GrossAnatomicalSite, pdf="figs/qcReport_IL003.pdf")
mdsPlot(RGset, n=1000, sampleGroups=targets$GrossAnatomicalSite)
mdsPlot(RGset, n=1000, sampleGroups=targets$GrossAnatomicalSite, sampleNames=sampleNames)



