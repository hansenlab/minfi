####################################################################################
## Create SNP annotation tables                                                   ##
## These tables can be cbinded to the IlluminaHumanMethylation450kanno data frame ##
## The rdas are in IlluminaHumanMethylation450kmanifest/data                      ##
####################################################################################
library(minfi)

## dbSNP 132
library(SNPlocs.Hsapiens.dbSNP.20101109)
snpAnno132 <- snpAnno()
# Add Illumina annotation SNPs since some appear to be missing from the SNPlocs package (e.g. rs34975136)
all(rownames(IlluminaHumanMethylation450kanno)==rownames(snpAnno132))
idx <- rowSums(IlluminaHumanMethylation450kanno[,c("Probe_SNPs", "Probe_SNPs_10")]!="")>0
snpAnno132$probeSnps[idx] <- TRUE
save(snpAnno132, file="snpAnno132.rda")

## dbSNP 131 
## v131 NOT YET ADDED TO THE MANIFEST PACKAGE. Didn't want to waste space since
## we'll probably reorganize soon
library(SNPlocs.Hsapiens.dbSNP.20100427)
snpAnno131 <- snpAnno()
save(snpAnno131, file="snpAnno131.rda")


