Benchmarking DelayedArray branch of **minfi**
================
Peter Hickey
11 April 2018

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
#>   2.764   0.421   3.223
system.time(MSetSWAN_v1.25.1 <- minfi::preprocessSWAN(RGSet))
#>    user  system elapsed 
#>  10.010   1.590  11.664
system.time(MSetIllumina_v1.25.1 <- minfi::detectionP(RGSet))
#>    user  system elapsed 
#>   4.629   0.568   5.220
system.time(MSetQuantile_v1.25.1 <- minfi::preprocessQuantile(RGSet))
#> [preprocessQuantile] Mapping to genome.
#> [preprocessQuantile] Fixing outliers.
#> [preprocessQuantile] Quantile normalizing.
#>    user  system elapsed 
#>  78.508   2.262  82.441
system.time(detP_v1.25.1 <- minfi::detectionP(RGSet))
#>    user  system elapsed 
#>   4.880   0.593   5.534

system.time(MSetFixedOutliers_v1.25.1 <- minfi::fixMethOutliers(MSet))
#>    user  system elapsed 
#>   0.666   0.203   0.872
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
```

Create the various *RGChannelSet* and *MethylSet* objects:

``` r
setRealizationBackend(NULL)
system.time(RGSet_in_memory_DelayedMatrix <- realize(RGSet))
#>    user  system elapsed 
#>   0.114   0.051   0.168
system.time(MSet_in_memory_DelayedMatrix <- realize(MSet))
#>    user  system elapsed 
#>   0.105   0.066   0.171

setRealizationBackend("HDF5Array")
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)
system.time(RGSet_on_disk_DelayedMatrix <- realize(RGSet))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.113   0.242   5.454
system.time(MSet_on_disk_DelayedMatrix <- realize(MSet))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.199   0.290   5.524

# Reset BACKEND and block size to defaults
setRealizationBackend(NULL)
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE)
```

### Using an ordinary *matrix*-backed *RGChannelSet*

``` r
system.time(MSetRaw <- preprocessRaw(RGSet))
#>    user  system elapsed 
#>   1.700   0.227   1.936
system.time(MSetIllumina <- preprocessIllumina(RGSet))
#>    user  system elapsed 
#>   2.840   0.563   3.415
system.time(MSetSWAN <- preprocessSWAN(RGSet))
#>    user  system elapsed 
#>   6.581   1.235   7.846
system.time(MSetQuantile <- preprocessQuantile(RGSet))
#> [preprocessQuantile] Mapping to genome.
#> [preprocessQuantile] Fixing outliers.
#> [preprocessQuantile] Quantile normalizing.
#>    user  system elapsed 
#>  72.603   1.480  74.929
system.time(detP <- detectionP(RGSet))
#>    user  system elapsed 
#>   4.944   0.654   5.729

system.time(MSetFixedOutliers <- fixMethOutliers(MSet))
#>    user  system elapsed 
#>   0.552   0.205   0.764
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
#>   6.113   0.956   7.225
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
#>  12.491   2.256  15.192
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
#>   7.446   1.459   9.063
system.time(
    MSetQuantile_in_memory_DelayedMatrix <- preprocessQuantile(
        RGSet_in_memory_DelayedMatrix))
#> [preprocessQuantile] Mapping to genome.
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
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
#>  83.032   3.684  88.210
system.time(
    detP_in_memory_DelayedMatrix <- detectionP(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   5.618   0.847   6.681

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
#>   1.457   0.751   2.523
```

#### Increased block size

``` r
setRealizationBackend(NULL)
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)

system.time(MSetRaw_in_memory_DelayedMatrix <- preprocessRaw(
    RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   1.652   0.221   1.912
system.time(
    MSetSWAN_in_memory_DelayedMatrix <- preprocessSWAN(
        RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   9.735   2.000  16.047
system.time(
    MSetIllumina_in_memory_DelayedMatrix <- preprocessIllumina(
        RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   3.430   0.701   5.094
system.time(
    MSetQuantile_in_memory_DelayedMatrix <- preprocessQuantile(
        RGSet_in_memory_DelayedMatrix))
#> [preprocessQuantile] Mapping to genome.
#> Processing block 1/1 ... OK
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
#>  79.604   2.887  88.954
system.time(
    detP_in_memory_DelayedMatrix <- detectionP(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.111   0.705   5.940

system.time(
    fixMethOutliers_in_memory_DelayedMatrix <- fixMethOutliers(
        MSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   0.801   0.298   1.148
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
#>  13.028   2.373  16.490
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
#>  22.283   4.076  27.319
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
#>  13.900   3.300  19.522
system.time(
    MSetQuantile_on_disk_DelayedMatrix <- preprocessQuantile(
        RGSet_on_disk_DelayedMatrix))
#> [preprocessQuantile] Mapping to genome.
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
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
#>  93.805   7.420 106.648
system.time(
    detP_on_disk_DelayedMatrix <- detectionP(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   8.516   1.985  11.239

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
#>   3.716   1.657   5.810
```

#### Increased block size

``` r
setRealizationBackend("HDF5Array")
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)

system.time(
    MSetRaw_on_disk_DelayedMatrix <- preprocessRaw(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   7.227   0.882  10.343
system.time(
    MSetSWAN_on_disk_DelayedMatrix <- preprocessSWAN(
        RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>  17.828   2.706  22.749
system.time(
    MSetIllumina_on_disk_DelayedMatrix <- preprocessIllumina(
        RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.562   1.403   7.051
system.time(
    MSetQuantile_on_disk_DelayedMatrix <- preprocessQuantile(
        RGSet_on_disk_DelayedMatrix))
#> [preprocessQuantile] Mapping to genome.
#> Processing block 1/1 ... OK
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
#>  83.292   3.850  91.163
system.time(
    detP_on_disk_DelayedMatrix <- detectionP(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.586   0.991   6.727

system.time(fixMethOutliers_in_memory_DelayedMatrix <-
                fixMethOutliers(MSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   1.471   0.575   2.070
```
