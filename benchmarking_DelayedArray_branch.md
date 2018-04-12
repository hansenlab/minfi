Benchmarking DelayedArray branch of **minfi**
================
Peter Hickey
12 April 2018

## Setup

``` r
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(minfiData))
```

## Using **minfi** v1.25.1

``` r
RGSet <- minfiData::RGsetEx
MSet <- minfiData::MsetEx

system.time(MSetRaw_v1.25.1 <- minfi::preprocessRaw(RGSet))
#>    user  system elapsed 
#>   2.645   0.401   3.067
system.time(MSetSWAN_v1.25.1 <- minfi::preprocessSWAN(RGSet))
#>    user  system elapsed 
#>  10.864   1.720  13.305
system.time(MSetIllumina_v1.25.1 <- minfi::detectionP(RGSet))
#>    user  system elapsed 
#>   4.518   0.560   5.109
system.time(MSetNoob_v1.25.1 <- minfi::preprocessNoob(RGSet))
#> [preprocessNoob] Applying R/G ratio flip to fix dye bias...
#>    user  system elapsed 
#> 100.740   3.325 108.946
system.time(detP_v1.25.1 <- minfi::detectionP(RGSet))
#>    user  system elapsed 
#>   5.008   0.665   5.799
system.time(oob_v1.25.1 <- minfi::getOOB(RGSet))
#>    user  system elapsed 
#>   0.740   0.115   0.873

system.time(MSetQuantile_v1.25.1 <- minfi::preprocessQuantile(MSet))
#> [preprocessQuantile] Mapping to genome.
#> [preprocessQuantile] Fixing outliers.
#> [preprocessQuantile] Quantile normalizing.
#>    user  system elapsed 
#>  76.424   1.988  84.967
system.time(MSetFixedOutliers_v1.25.1 <- minfi::fixMethOutliers(MSet))
#>    user  system elapsed 
#>   0.532   0.141   0.680
```

## Using DelayedArray branch of **minfi**

Some setting up:

``` r
source("setup_devel_workspace.R")
#> Loading required package: rhdf5
#> 
#> Attaching package: 'DelayedMatrixStats'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyMissings, colAnyNAs, colAnys, colAvgsPerRowSet,
#>     colCollapse, colCounts, colCummaxs, colCummins, colCumprods,
#>     colCumsums, colDiffs, colIQRDiffs, colIQRs, colLogSumExps,
#>     colMadDiffs, colMads, colMeans2, colMedians, colOrderStats,
#>     colProds, colQuantiles, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyMissings, rowAnyNAs, rowAnys,
#>     rowAvgsPerColSet, rowCollapse, rowCounts, rowCummaxs,
#>     rowCummins, rowCumprods, rowCumsums, rowDiffs, rowIQRDiffs,
#>     rowIQRs, rowLogSumExps, rowMadDiffs, rowMads, rowMeans2,
#>     rowMedians, rowOrderStats, rowProds, rowQuantiles, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs,
#>     rowVars, rowWeightedMads, rowWeightedMeans,
#>     rowWeightedMedians, rowWeightedSds, rowWeightedVars
#> The following object is masked from 'package:Biobase':
#> 
#>     rowMedians
#> Warning in .removePreviousCoerce(class1, class2, where, prevIs): methods
#> currently exist for coercing from "MethylSet" to "Vector"; they will be
#> replaced.
```

Create the various *RGChannelSet* and *MethylSet* objects:

``` r
setRealizationBackend(NULL)
system.time(RGSet_in_memory_DelayedMatrix <- realize(RGSet))
#> Found more than one class "RGChannelSet" in cache; using the first, from namespace 'minfi'
#> Also defined by '.GlobalEnv'
#> Found more than one class "RGChannelSet" in cache; using the first, from namespace 'minfi'
#> Also defined by '.GlobalEnv'
#> Found more than one class "RGChannelSet" in cache; using the first, from namespace 'minfi'
#> Also defined by '.GlobalEnv'
#> Found more than one class "RGChannelSet" in cache; using the first, from namespace 'minfi'
#> Also defined by '.GlobalEnv'
#> Found more than one class "RGChannelSet" in cache; using the first, from namespace 'minfi'
#> Also defined by '.GlobalEnv'
#>    user  system elapsed 
#>   0.109   0.046   0.157
system.time(MSet_in_memory_DelayedMatrix <- realize(MSet))
#> Found more than one class "MethylSet" in cache; using the first, from namespace 'minfi'
#> Also defined by '.GlobalEnv'
#> Found more than one class "MethylSet" in cache; using the first, from namespace 'minfi'
#> Also defined by '.GlobalEnv'
#> Found more than one class "MethylSet" in cache; using the first, from namespace 'minfi'
#> Also defined by '.GlobalEnv'
#> Found more than one class "MethylSet" in cache; using the first, from namespace 'minfi'
#> Also defined by '.GlobalEnv'
#>    user  system elapsed 
#>   0.146   0.098   0.258

setRealizationBackend("HDF5Array")
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)
system.time(RGSet_on_disk_DelayedMatrix <- realize(RGSet))
#> Processing block 1/1 ...
#> OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.156   0.279   5.547
system.time(MSet_on_disk_DelayedMatrix <- realize(MSet))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.249   0.254   5.677

# Reset BACKEND and block size to defaults
setRealizationBackend(NULL)
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE)
```

### Using an ordinary *matrix*-backed *RGChannelSet*

``` r
system.time(MSetRaw <- preprocessRaw(RGSet))
#>    user  system elapsed 
#>   1.810   0.233   2.144
system.time(MSetIllumina <- preprocessIllumina(RGSet))
#> Found more than one class "RGChannelSet" in cache; using the first, from namespace 'minfi'
#> Also defined by '.GlobalEnv'
#>    user  system elapsed 
#>   3.175   0.731   4.061
system.time(MSetSWAN <- preprocessSWAN(RGSet))
#>    user  system elapsed 
#>   7.193   1.370   8.730
system.time(MSetNoob <- preprocessNoob(RGSet))
#> [dyeCorrection] Applying R/G ratio flip to fix dye bias
#>    user  system elapsed 
#> 100.006   3.236 109.873
system.time(detP <- detectionP(RGSet))
#>    user  system elapsed 
#>   4.560   0.573   5.164
system.time(oob <- getOOB(RGSet))
#>    user  system elapsed 
#>   0.592   0.074   0.668

system.time(MSetQuantile <- preprocessQuantile(MSet))
#> [preprocessQuantile] Mapping to genome.
#> Found more than one class "MethylSet" in cache; using the first, from namespace 'minfi'
#> Also defined by '.GlobalEnv'
#> [preprocessQuantile] Fixing outliers.
#> [preprocessQuantile] Quantile normalizing.
#>    user  system elapsed 
#>  73.024   1.514  78.194
system.time(MSetFixedOutliers <- fixMethOutliers(MSet))
#>    user  system elapsed 
#>   0.520   0.170   0.691
```

### Using an in-memory *DelayedMatrix*-backed *RGChannelSet*

#### Default block size

``` r
setRealizationBackend(NULL)
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE)

system.time(
    MSetRaw_in_memory_DelayedMatrix <- preprocessRaw(
        RGSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   5.317   0.783   6.115
system.time(
    MSetSWAN_in_memory_DelayedMatrix <- preprocessSWAN(
        RGSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/3 ... OK
#> Processing block 2/3 ... OK
#> Processing block 3/3 ... OK
#> Processing block 1/3 ... OK
#> Processing block 2/3 ... OK
#> Processing block 3/3 ... OK
#>    user  system elapsed 
#>  11.558   2.031  13.667
system.time(
    MSetIllumina_in_memory_DelayedMatrix <- preprocessIllumina(
        RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   6.985   1.334   8.397
system.time(
    MSetNoob_in_memory_DelayedMatrix <- preprocessNoob(
        RGSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/3 ... Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> [dyeCorrection] Applying R/G ratio flip to fix dye bias
#> OK
#> Processing block 2/3 ... Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> [dyeCorrection] Applying R/G ratio flip to fix dye bias
#> OK
#> Processing block 3/3 ... Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> [dyeCorrection] Applying R/G ratio flip to fix dye bias
#> OK
#>    user  system elapsed 
#> 100.128   4.446 105.707
system.time(
    detP_in_memory_DelayedMatrix <- detectionP(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   6.074   0.937   7.422
system.time(
    oob_in_memory_DelayedMatrix <- getOOB(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   2.553   0.364   3.607

system.time(
    MSetQuantile_in_memory_DelayedMatrix <- preprocessQuantile(
        MSet_in_memory_DelayedMatrix))
#> [preprocessQuantile] Mapping to genome.
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> [preprocessQuantile] Fixing outliers.
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> [preprocessQuantile] Quantile normalizing.
#>    user  system elapsed 
#>  72.975   2.379  78.508
system.time(
    fixMethOutliers_in_memory_DelayedMatrix <- fixMethOutliers(
        MSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   1.434   0.656   2.113
```

#### Increased block size

``` r
setRealizationBackend(NULL)
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)

system.time(MSetRaw_in_memory_DelayedMatrix <- preprocessRaw(
    RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   1.570   0.211   1.800
system.time(
    MSetSWAN_in_memory_DelayedMatrix <- preprocessSWAN(
        RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   7.524   1.440   9.794
system.time(
    MSetIllumina_in_memory_DelayedMatrix <- preprocessIllumina(
        RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   2.655   0.527   3.266
system.time(
    MSetNoob_in_memory_DelayedMatrix <- preprocessNoob(
        RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> [dyeCorrection] Applying R/G ratio flip to fix dye bias
#> OK
#>    user  system elapsed 
#>  98.265   3.477 106.131
system.time(
    detP_in_memory_DelayedMatrix <- detectionP(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   4.986   0.769   5.824
system.time(
    oob_in_memory_DelayedMatrix <- getOOB(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   0.841   0.122   1.014

system.time(
    MSetQuantile_in_memory_DelayedMatrix <- preprocessQuantile(
        MSet_in_memory_DelayedMatrix))
#> [preprocessQuantile] Mapping to genome.
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> [preprocessQuantile] Fixing outliers.
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> [preprocessQuantile] Quantile normalizing.
#>    user  system elapsed 
#>  71.589   1.887  75.841
system.time(
    fixMethOutliers_in_memory_DelayedMatrix <- fixMethOutliers(
        MSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   0.762   0.258   1.024
```

### Using an on-disk (HDF5) *DelayedMatrix*-backed *RGChannelSet*

Using the **HDF5Array** backend

#### Default block size

``` r
setRealizationBackend("HDF5Array")
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE)

system.time(
    MSetRaw_on_disk_DelayedMatrix <- preprocessRaw(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>  11.523   2.136  13.808
system.time(
    MSetSWAN_on_disk_DelayedMatrix <- preprocessSWAN(
        RGSet_on_disk_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/3 ... OK
#> Processing block 2/3 ... OK
#> Processing block 3/3 ... OK
#> Processing block 1/3 ... OK
#> Processing block 2/3 ... OK
#> Processing block 3/3 ... OK
#>    user  system elapsed 
#>  21.347   3.813  25.418
system.time(
    MSetIllumina_on_disk_DelayedMatrix <- preprocessIllumina(
        RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>  11.873   3.003  15.252
system.time(
    MSetNoob_on_disk_DelayedMatrix <- preprocessNoob(
        RGSet_on_disk_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/3 ... Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> [dyeCorrection] Applying R/G ratio flip to fix dye bias
#> OK
#> Processing block 2/3 ... Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> [dyeCorrection] Applying R/G ratio flip to fix dye bias
#> OK
#> Processing block 3/3 ... Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> [dyeCorrection] Applying R/G ratio flip to fix dye bias
#> OK
#>    user  system elapsed 
#> 118.339   8.181 130.714
system.time(
    detP_on_disk_DelayedMatrix <- detectionP(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   9.446   2.133  13.221
system.time(oob_on_disk_DelayedMatrix <- getOOB(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   5.006   0.506   5.766

system.time(
    MSetQuantile_on_disk_DelayedMatrix <- preprocessQuantile(
        MSet_on_disk_DelayedMatrix))
#> [preprocessQuantile] Mapping to genome.
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> [preprocessQuantile] Fixing outliers.
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> [preprocessQuantile] Quantile normalizing.
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>  82.523   5.652  91.085
system.time(
    fixMethOutliers_in_memory_DelayedMatrix <- fixMethOutliers(
        MSet_on_disk_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   3.321   1.480   4.843
```

#### Increased block size

``` r
setRealizationBackend("HDF5Array")
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)

system.time(
    MSetRaw_on_disk_DelayedMatrix <- preprocessRaw(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.776   0.563   6.441
system.time(
    MSetSWAN_on_disk_DelayedMatrix <- preprocessSWAN(
        RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>  15.772   2.361  18.436
system.time(
    MSetIllumina_on_disk_DelayedMatrix <- preprocessIllumina(
        RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   7.397   1.697  11.466
system.time(
    MSetNoob_on_disk_DelayedMatrix <- preprocessNoob(
        RGSet_on_disk_DelayedMatrix))
#> Warning in h5createDataset(filepath, name, dim, storage.mode = type, size
#> = size, : You created a large dataset with compression and chunking. The
#> chunk size is equal to the dataset dimensions. If you want to read subsets
#> of the dataset, you should test smaller chunk sizes to improve read times.
#> Turn off this warning with showWarnings=FALSE.
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Visiting block 1/1 ... OK
#> [dyeCorrection] Applying R/G ratio flip to fix dye bias
#> OK
#>    user  system elapsed 
#> 108.790   5.116 120.332
system.time(
    detP_on_disk_DelayedMatrix <- detectionP(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.763   0.972   6.933
system.time(oob_on_disk_DelayedMatrix <- getOOB(RGSet_in_memory_DelayedMatrix))
#> Warning in h5createDataset(filepath, name, dim, storage.mode = type, size
#> = size, : You created a large dataset with compression and chunking. The
#> chunk size is equal to the dataset dimensions. If you want to read subsets
#> of the dataset, you should test smaller chunk sizes to improve read times.
#> Turn off this warning with showWarnings=FALSE.
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   2.283   0.182   2.507

system.time(
    MSetQuantile_on_disk_DelayedMatrix <- preprocessQuantile(
        MSet_on_disk_DelayedMatrix))
#> [preprocessQuantile] Mapping to genome.
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> [preprocessQuantile] Fixing outliers.
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> [preprocessQuantile] Quantile normalizing.
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>  77.520   3.481  83.414
system.time(fixMethOutliers_in_memory_DelayedMatrix <-
                fixMethOutliers(MSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   1.496   0.506   2.253
```
