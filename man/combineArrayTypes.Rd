\name{combineArrayTypes}
\alias{combineArrayTypes}
\title{
A method for combining different types of methylation arrays into a
virtual array.
}
\description{
A method for combining different types of methylation arrays into a
virtual array. Specifically, the 450k array and the EPIC array share
many probes in common.  This functions combines data from the two
different array types and outputs a data object of the user specified
type.  Essentially this new object will be like (for example) an EPIC
array with many probes missing.
}
\usage{
combineArrayTypes(rgSet1, rgSet2, outType = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC"), verbose = TRUE)
}
\arguments{
  \item{rgSet1}{The first \code{RGChannelSet}.}
  \item{rgSet2}{The second \code{RGChannelSet}.}
  \item{outType}{The array type of the output.}
  \item{verbose}{Should the function be verbose?}
}
\details{
  Still under development.

  SNP probes. Control probes.
}
\value{
  An object of class \code{RGChannelSet}, with the type of array given
  by the \code{outType} argument.
}
\author{
  Jean-Philippe Fortin and Kasper D. Hansen.
}
\examples{
if(require(minfiData) && require(minfiDataEPIC)) {
    library(minfiData)
    library(minfiDataEPIC)
    data(RGsetEx)
    data(RSsetEPIC)
    rgSet <- combineArrayTypes(RGsetEPIC, RGsetEx)
    rgSet
}
}
