Benchmarking DelayedArray branch of **minfi**
================
Peter Hickey
3 April 2018

## Setup

``` r
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(minfiData))
```

## Using **minfi** v1.25.1

``` r
RGSet <- minfiData::RGsetEx
system.time(MSetRaw_v1.25.1 <- minfi::preprocessRaw(RGSet))
#>    user  system elapsed 
#>   3.594   0.522   4.377
system.time(MSetIllumina_v1.25.1 <- minfi::detectionP(RGSet))
#>    user  system elapsed 
#>   6.170   0.600   6.867
system.time(detP_v1.25.1 <- minfi::detectionP(RGSet))
#>    user  system elapsed 
#>   5.841   0.708   6.831

MSet <- minfiData::MsetEx
system.time(MSetFixedOutliers_v1.25.1 <- minfi::fixMethOutliers(MSet))
#>    user  system elapsed 
#>   0.698   0.254   0.992
```

## Using DelayedArray branch of **minfi**

Some setting up:

``` r
suppressPackageStartupMessages(library(DelayedArray))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(DelayedMatrixStats))

source("R/utils.R")
source("R/DelayedArray_utils.R")
source("R/preprocessRaw.R")
source("R/preprocessIllumina.R")
source("R/detectionP.R")
source("R/minfiQC.R")

DelayedArray:::set_verbose_block_processing(TRUE)
#> [1] FALSE

DEFAULT_BLOCK_SIZE <- getOption("DelayedArray.block.size")
DEFAULT_BLOCK_SIZE
#> [1] 4500000
```

Create the various *RGChannelSet* and *MethylSet* objects:

``` r
RGSet_in_memory_DelayedMatrix <- minfiData::RGsetEx
system.time(
    assays(RGSet_in_memory_DelayedMatrix) <- endoapply(
        assays(RGSet_in_memory_DelayedMatrix),
        DelayedArray)
)
#>    user  system elapsed 
#>   0.032   0.009   0.042

RGSet_on_disk_DelayedMatrix <- minfiData::RGsetEx
system.time(
    assays(RGSet_on_disk_DelayedMatrix) <- endoapply(
        assays(RGSet_on_disk_DelayedMatrix),
        writeHDF5Array))
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
#>   7.184   1.374   8.943

MSet_in_memory_DelayedMatrix <- minfiData::MsetEx
system.time(
    assays(MSet_in_memory_DelayedMatrix) <- endoapply(
        assays(MSet_in_memory_DelayedMatrix),
        DelayedArray)
)
#>    user  system elapsed 
#>   0.026   0.017   0.045

MSet_on_disk_DelayedMatrix <- minfiData::MsetEx
system.time(
    assays(MSet_on_disk_DelayedMatrix) <- endoapply(
        assays(MSet_on_disk_DelayedMatrix),
        writeHDF5Array))
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
#>   6.038   1.260   7.537
```

### Using an ordinary *matrix*-backed *RGChannelSet*

``` r
system.time(MSetRaw <- preprocessRaw(RGSet))
#>    user  system elapsed 
#>   1.957   0.256   2.246
system.time(MSetIllumina <- preprocessIllumina(RGSet))
#>    user  system elapsed 
#>   4.603   0.792   5.865
system.time(detP <- detectionP(RGSet))
#>    user  system elapsed 
#>   6.511   0.913   8.281

system.time(MSetFixedOutliers <- fixMethOutliers(MSet))
#>    user  system elapsed 
#>   0.644   0.254   0.947
```

### Using an in-memory *DelayedMatrix*-backed *RGChannelSet*

#### Default block size

``` r
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE)
system.time(
    MSetRaw_in_memory_DelayedMatrix <- preprocessRaw(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   7.305   1.078   8.880
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
#>   9.173   1.640  11.431
system.time(
    detP_in_memory_DelayedMatrix <- detectionP(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   6.422   0.974   7.615

system.time(fixMethOutliers_in_memory_DelayedMatrix <-
                fixMethOutliers(MSet_in_memory_DelayedMatrix))
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
#>   1.290   0.641   1.979
```

#### Increased block size

``` r
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)
system.time(MSetRaw_in_memory_DelayedMatrix <- preprocessRaw(
    RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   2.035   0.333   2.487
system.time(MSetIllumina_in_memory_DelayedMatrix <- preprocessIllumina(
    RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   3.239   0.706   4.088
system.time(
    detP_in_memory_DelayedMatrix <- detectionP(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.877   0.749   6.937

system.time(fixMethOutliers_in_memory_DelayedMatrix <-
                fixMethOutliers(MSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   1.019   0.403   1.500
```

### Using an on-disk (HDF5) *DelayedMatrix*-backed *RGChannelSet*

Using the **HDF5Array** backend:

``` r
setRealizationBackend("HDF5Array")
```

#### Default block size

``` r
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE)
system.time(
    MSet_on_disk_DelayedMatrix <- preprocessRaw(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>  13.990   2.707  18.463
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
#>  15.043   3.825  25.527
system.time(
   detP_on_disk_DelayedMatrix <- detectionP(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   8.400   1.791  10.585

system.time(fixMethOutliers_in_memory_DelayedMatrix <-
                fixMethOutliers(MSet_on_disk_DelayedMatrix))
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
#>   3.819   1.708   6.051
```

#### Increased block size

``` r
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)
system.time(
    MSet_on_disk_DelayedMatrix <- preprocessRaw(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   6.975   0.783   9.788
system.time(
    MSetIllumina_on_disk_DelayedMatrix <- preprocessIllumina(
        RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   6.900   1.618   9.888
system.time(
   detP_on_disk_DelayedMatrix <- detectionP(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   6.729   1.079   8.955

system.time(fixMethOutliers_in_memory_DelayedMatrix <-
                fixMethOutliers(MSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   1.535   0.558   2.122
```
