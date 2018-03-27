Benchmarking DelayedArray branch of **minfi**
================
Peter Hickey
27 March 2018

## Setup

``` r
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(minfiData))
```

## Using **minfi** v1.25.1

``` r
RGSet <- minfiData::RGsetEx
system.time(MSetRaw_v1.25.1 <- preprocessRaw(RGSet))
#>    user  system elapsed 
#>   3.540   0.506   4.102
system.time(MSetIllumina_v1.25.1 <- detectionP(RGSet))
#>    user  system elapsed 
#>   6.374   0.676   7.167
system.time(detP_v1.25.1 <- detectionP(RGSet))
#>    user  system elapsed 
#>   7.072   0.723   7.966
```

## Using DelayedArray branch of **minfi**

Some setting up:

``` r
suppressPackageStartupMessages(library(DelayedArray))
suppressPackageStartupMessages(library(HDF5Array))

source("R/utils.R")
source("R/DelayedArray_utils.R")
source("R/preprocessRaw.R")
source("R/preprocessIllumina.R")
source("R/detectionP.R")

DelayedArray:::set_verbose_block_processing(TRUE)
#> [1] FALSE

DEFAULT_BLOCK_SIZE <- getOption("DelayedArray.block.size")
DEFAULT_BLOCK_SIZE
#> [1] 4500000
```

Create the various *RGChannelSet* objects:

``` r
RGSet_in_memory_DelayedMatrix <- minfiData::RGsetEx
system.time(
    assays(RGSet_in_memory_DelayedMatrix) <- endoapply(
        assays(RGSet_in_memory_DelayedMatrix),
        DelayedArray)
)
#>    user  system elapsed 
#>   0.050   0.016   0.067

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
#>   7.980   0.982   9.170
```

### Using an ordinary *matrix*-backed *RGChannelSet*

``` r
system.time(MSetRaw <- preprocessRaw(RGSet))
#>    user  system elapsed 
#>   2.139   0.279   2.498
system.time(MSetIllumina <- preprocessIllumina(RGSet))
#>    user  system elapsed 
#>   4.286   0.708   5.140
system.time(detP <- detectionP(RGSet))
#>    user  system elapsed 
#>   5.869   0.756   6.771
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
#>   7.579   1.082   8.905
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
#>  11.177   1.992  13.858
system.time(
    detP_in_memory_DelayedMatrix <- detectionP(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   8.840   1.380  16.223
```

#### Increased block size

``` r
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)
system.time(MSetRaw_in_memory_DelayedMatrix <- preprocessRaw(
    RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   3.229   0.459   6.517
system.time(MSetIllumina_in_memory_DelayedMatrix <- preprocessIllumina(
    RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   5.238   0.941  10.055
system.time(
    detP_in_memory_DelayedMatrix <- detectionP(RGSet_in_memory_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   7.085   1.032  12.026
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
#>  14.006   2.481  18.932
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
#>  15.236   3.342  19.883
system.time(
   detP_on_disk_DelayedMatrix <- detectionP(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/6 ... OK
#> Processing block 2/6 ... OK
#> Processing block 3/6 ... OK
#> Processing block 4/6 ... OK
#> Processing block 5/6 ... OK
#> Processing block 6/6 ... OK
#>    user  system elapsed 
#>   9.897   2.157  13.422
```

#### Increased block size

``` r
options(DelayedArray.block.size = DEFAULT_BLOCK_SIZE * 100L)
system.time(
    MSet_on_disk_DelayedMatrix <- preprocessRaw(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   6.674   0.627   7.733
system.time(
    MSetIllumina_on_disk_DelayedMatrix <- preprocessIllumina(
        RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   8.280   1.429  12.522
system.time(
   detP_on_disk_DelayedMatrix <- detectionP(RGSet_on_disk_DelayedMatrix))
#> Processing block 1/1 ... OK
#>    user  system elapsed 
#>   6.722   1.056   8.199
```
