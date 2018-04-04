Benchmarking DelayedArray branch of **minfi**
================
Peter Hickey
4 April 2018

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
#>   3.666   0.534   4.922
system.time(MSetSWAN_v1.25.1 <- minfi::preprocessSWAN(RGSet))
#>    user  system elapsed 
#>  12.472   1.998  16.564
system.time(MSetIllumina_v1.25.1 <- minfi::detectionP(RGSet))
#>    user  system elapsed 
#>   5.275   0.708   6.354
system.time(detP_v1.25.1 <- minfi::detectionP(RGSet))
#>    user  system elapsed 
#>   5.618   0.741   7.660

MSet <- minfiData::MsetEx
system.time(MSetFixedOutliers_v1.25.1 <- minfi::fixMethOutliers(MSet))
#>    user  system elapsed 
#>   0.645   0.246   0.975
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
source("R/preprocessSwan.R")

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
#>   0.024   0.008   0.032

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
#>   6.811   1.209   8.624

MSet_in_memory_DelayedMatrix <- minfiData::MsetEx
system.time(
    assays(MSet_in_memory_DelayedMatrix) <- endoapply(
        assays(MSet_in_memory_DelayedMatrix),
        DelayedArray)
)
#>    user  system elapsed 
#>   0.025   0.019   0.044

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
#>   5.826   1.187   7.177
```

### Using an ordinary *matrix*-backed *RGChannelSet*

``` r
system.time(MSetRaw <- preprocessRaw(RGSet))
#>    user  system elapsed 
#>   1.683   0.225   1.925
system.time(MSetIllumina <- preprocessIllumina(RGSet))
#>    user  system elapsed 
#>   3.388   0.608   4.089
system.time(MSetSWAN <- preprocessSWAN(RGSet))
#>    user  system elapsed 
#>   6.713   1.331   8.190
system.time(detP <- detectionP(RGSet))
#>    user  system elapsed 
#>   5.612   0.777   6.790

system.time(MSetFixedOutliers <- fixMethOutliers(MSet))
#>    user  system elapsed 
#>   0.656   0.249   0.926
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
#>   8.326   1.277  12.519
system.time(MSetSWAN_in_memory_DelayedMatrix <- preprocessSWAN(
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
#>    user  system elapsed 
#>  15.022   2.809  20.762
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
#>   7.355   1.523   9.078
system.time(
    detP_in_memory_DelayedMatrix <- detectionP(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   5.790   0.912   6.788

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
#>    user  system elapsed 
#>   1.157   0.561   1.736
```

#### Increased block size

``` r
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)
system.time(MSetRaw_in_memory_DelayedMatrix <- preprocessRaw(
    RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   1.745   0.264   2.033
system.time(
    MSetSWAN_in_memory_DelayedMatrix <- preprocessSWAN(
        RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   6.878   1.347   8.335
system.time(
    MSetIllumina_in_memory_DelayedMatrix <- preprocessIllumina(
        RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   3.472   0.568   4.097
system.time(
    detP_in_memory_DelayedMatrix <- detectionP(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.165   0.702   6.233

system.time(
    fixMethOutliers_in_memory_DelayedMatrix <- fixMethOutliers(
        MSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   0.990   0.557   1.908
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
    MSetRaw_on_disk_DelayedMatrix <- preprocessRaw(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>  11.590   2.575  14.662
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
#>  25.274   5.875  33.219

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
#>  12.023   3.219  15.729
system.time(
    detP_on_disk_DelayedMatrix <- detectionP(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   8.463   2.002  11.110

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
#>   4.154   2.051   6.653
```

#### Increased block size

``` r
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)
system.time(
    MSetRaw_on_disk_DelayedMatrix <- preprocessRaw(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   6.227   0.787   7.312
system.time(
    MSetSWAN_on_disk_DelayedMatrix <- preprocessSWAN(
        RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>  16.869   2.800  20.560
system.time(
    MSetIllumina_on_disk_DelayedMatrix <- preprocessIllumina(
        RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   6.795   1.589  10.046
system.time(
    detP_on_disk_DelayedMatrix <- detectionP(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   7.401   1.255  10.741

system.time(fixMethOutliers_in_memory_DelayedMatrix <-
                fixMethOutliers(MSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   1.550   0.511   2.079
```
