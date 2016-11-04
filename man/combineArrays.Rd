\name{combineArrays}
\alias{combineArrays}
\alias{combineArrays,RGChannelSet,RGChannelSet-method}
\alias{combineArrays,MethylSet,MethylSet-method}
\alias{combineArrays,GenomicMethylSet,GenomicMethylSet-method}
\alias{combineArrays,RatioSet,RatioSet-method}
\alias{combineArrays,GenomicRatioSet,GenomicRatioSet-method}
\title{
A method for combining different types of methylation arrays into a
virtual array.
}
\description{
  A method for combining different types of methylation arrays into a
virtual array. The three generations of Illumina methylation arrays are
supported: the 27k, the 450k and the EPIC arrays.  Specifically, the
450k array and the EPIC array share many probes in common.  This
function combines data from the two different array types and outputs a
data object of the user-specified type.  Essentially, this new object
will be like (for example) an EPIC array with many probes missing.
}
\usage{
\S4method{combineArrays}{RGChannelSet,RGChannelSet}(object1, object2,
                  outType = c("IlluminaHumanMethylation450k",
                              "IlluminaHumanMethylationEPIC"),
                  verbose = TRUE)
\S4method{combineArrays}{MethylSet,MethylSet}(object1, object2,
                  outType = c("IlluminaHumanMethylation450k",
                              "IlluminaHumanMethylationEPIC",
                              "IlluminaHumanMethylation27k"),
                  verbose = TRUE)
\S4method{combineArrays}{RatioSet,RatioSet}(object1, object2,
                  outType = c("IlluminaHumanMethylation450k",
                              "IlluminaHumanMethylationEPIC",
                              "IlluminaHumanMethylation27k"),
                  verbose = TRUE)
\S4method{combineArrays}{GenomicMethylSet,GenomicMethylSet}(object1, object2,
                  outType = c("IlluminaHumanMethylation450k",
                              "IlluminaHumanMethylationEPIC",
                              "IlluminaHumanMethylation27k"),
                  verbose = TRUE)
\S4method{combineArrays}{GenomicRatioSet,GenomicRatioSet}(object1, object2,
                  outType = c("IlluminaHumanMethylation450k",
                              "IlluminaHumanMethylationEPIC",
                              "IlluminaHumanMethylation27k"),
                  verbose = TRUE)
}
\arguments{
  \item{object1}{The first object.}
  \item{object2}{The second object.}
  \item{outType}{The array type of the output.}
  \item{verbose}{Should the function be verbose?}
}
\details{
  FIXME: describe the RCChannelSet combination.
}
\value{
  The output object has the same class as the two input objects, that is
  either an \code{RGChannelSet}, a \code{MethylSet}, a \code{RatioSet}, a
  \code{GemomicMethylSet} or a \code{GenomicRatioSet}, with the type of
  the array given by the \code{outType} argument.
}
\author{
  Jean-Philippe Fortin and Kasper D. Hansen.
}
\examples{
if(require(minfiData) && require(minfiDataEPIC)) {
  data(RGsetEx.sub)
  data(RGsetEPIC)
  rgSet <- combineArrays(RGsetEPIC, RGsetEx.sub)
  rgSet
}
}
