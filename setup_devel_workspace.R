# Load everything needed to start devel work in the DelayedArray branch of
# minfi. Can't simply use `devtools::load_all()` because of some sort of
# circular dependency amongst minfi's dependencies.
library(minfiData)
library(HDF5Array)
library(profvis)
library(profmem)
library(DelayedMatrixStats)

DelayedArray:::set_verbose_block_processing(TRUE)

# Add files here as they've been updated to handle DelayedArray objects
source("R/DelayedArray_utils.R")
source("R/utils.R")
source("R/preprocessRaw.R")
source("R/preprocessIllumina.R")
source("R/detectionP.R")
source("R/qc.R")
source("R/read.meth.R")
source("R/getSex.R")
source("R/minfiQC.R")
source("R/preprocessQuantile.R")
source("R/preprocessSwan.R")
