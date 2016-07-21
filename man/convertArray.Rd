\name{convertArray}
\alias{convertArray}
\alias{convertArray,RGChannelSet-method}
\alias{convertArray,MethylSet-method}
\alias{convertArray,GenomicMethylSet-method}
\alias{convertArray,RatioSet-method}
\alias{convertArray,GenomicRatioSet-method}
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
\S4method{convertArray}{RGChannelSet}(object,
                  outType = c("IlluminaHumanMethylation450k",
                              "IlluminaHumanMethylationEPIC"),
                  verbose = TRUE)
\S4method{convertArray}{MethylSet}(object,
                  outType = c("IlluminaHumanMethylation450k",
                              "IlluminaHumanMethylationEPIC",
                              "IlluminaHumanMethylation27k"),
                  verbose = TRUE)
\S4method{convertArray}{RatioSet}(object,
                  outType = c("IlluminaHumanMethylation450k",
                              "IlluminaHumanMethylationEPIC",
                              "IlluminaHumanMethylation27k"),
                  verbose = TRUE)
\S4method{convertArray}{GenomicMethylSet}(object, 
                  outType = c("IlluminaHumanMethylation450k",
                              "IlluminaHumanMethylationEPIC",
                              "IlluminaHumanMethylation27k"),
                  verbose = TRUE)
\S4method{convertArray}{GenomicRatioSet}(object, 
                  outType = c("IlluminaHumanMethylation450k",
                              "IlluminaHumanMethylationEPIC",
                              "IlluminaHumanMethylation27k"),
                  verbose = TRUE)
}
\arguments{
  \item{object}{The input object.}
  \item{outType}{The array type of the output.}
  \item{verbose}{Should the function be verbose?}
}
\details{
  Still under development.

  SNP probes. Control probes.
}
\value{
  FIXME
}
\author{
  Jean-Philippe Fortin and Kasper D. Hansen.
}
\examples{
if(require(minfiData)) {
    library(minfiData)
    data(RGsetEx)
    rgSet <- convertArray(RGsetEx, outType = "IlluminaHumanMethylationEPIC")
    rgSet
}
}
