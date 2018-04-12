# Load everything needed to start devel work in the DelayedArray branch of
# minfi. Can't simply use `devtools::load_all()` because of some sort of
# circular dependency amongst minfi's dependencies.
suppressPackageStartupMessages(library(minfiData))
library(HDF5Array)
library(profvis)
library(profmem)
library(DelayedMatrixStats)

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
source("R/RGChannelSet-class.R")
source("R/RGChannelSetExtended-class.R")
source("R/preprocessNoob.R")
source("R/MethylSet-class.R")

# These functions are properly imported by the package but not visible when the
# workspace is set up using this script
normalize.quantiles <- preprocessCore::normalize.quantiles
normalize.quantiles.use.target <- preprocessCore::normalize.quantiles.use.target
type <- DelayedArray::type
huber <- MASS::huber
normexp.signal <- limma::normexp.signal

# These methods are properly defined by the package but not visible when the
# workspace is set up using this script
setMethod("mapToGenome", signature(object = "RGChannelSet"),
          function(object, ...) {
              object <- preprocessRaw(object)
              callGeneric(object, ...)
          })

# Some DelayedArray functionality
DelayedArray:::set_verbose_block_processing(TRUE)
DEFAULT_BLOCK_SIZE <- getOption("DelayedArray.block.size")
DEFAULT_BLOCK_SIZE
