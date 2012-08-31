## setClass("MethylSetGenome",
##          contains = "SummarizedExperiment")

## library(minfiLocal)
## data(MsetEx)


## sset <- SummarizedExperiment(assays = SimpleList(Meth = getMeth(MsetEx),
##                              Unmeth = getUnmeth(MsetEx)),
##                              colData = phenoData(MsetEx),
##                              rowData = getLocations(MsetEx, genomeBuild = "hg19", returnAs = "GRanges"))
                             
