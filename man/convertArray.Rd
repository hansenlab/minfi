\name{convertArray}
\alias{convertArray}
\alias{convertArray,RGChannelSet-method}
\alias{convertArray,MethylSet-method}
\alias{convertArray,GenomicMethylSet-method}
\alias{convertArray,RatioSet-method}
\alias{convertArray,GenomicRatioSet-method}
\title{
A method for converting a type of methylation arrays into a
virtual array of another type.
}
\description{
  A method for converting a type of methylation array into a array of
another type. The three generations of Illumina methylation arrays are
supported: the 27k, the 450k and the EPIC arrays.  Specifically, the
450k array and the EPIC array share many probes in common. For
\code{RGChannelSet}, this function will convert an EPIC array into a
450k array (or vice-versa) by dropping probes that differ between the
two arrays. Because most of the probes on the 27k array have a different
chemistry than the 450k and EPIC probes, converting an 27k
\code{RGChannelSet} into another array is not supported.  Each array can
be converted into another array at the CpG site level, that is any
\code{MethylSet} and \code{RatioSet} (or \code{GenomicMethylSet} and
\code{GenomicRatioSet}) can be converted to a 27k, 450k or EPIC array.
The output array is specified by the \code{outType} argument.
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
  FIXME: describe the RGChannelSet conversion.
}
\value{
  The output object has the same class as the input object, that is
  either an \code{RGChannelSet}, a \code{MethylSet}, a \code{RatioSet}, a
  \code{GemomicMethylSet} or a \code{GenomicRatioSet}, with the type of
  the array given by the \code{outType} argument.
  }
\author{
  Jean-Philippe Fortin and Kasper D. Hansen.
}
\examples{
if(require(minfiData)) {
  data(RGsetEx.sub)
  rgSet <- convertArray(RGsetEx.sub, outType = "IlluminaHumanMethylationEPIC")
  rgSet
}
}
