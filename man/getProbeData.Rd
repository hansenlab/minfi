\name{getProbeData}
\alias{getProbeData}
\alias{getProbeInfo}
\alias{getManifestInfo}
\alias{getControlAddress}
\title{
  Utility functions for retrieving array design information for Illumina
  methylation microarrays.
}
\description{
  A set of functions for retrieving array design information for
  Illumina methylation microarrays.
}
\usage{
getProbeData(object)
getProbeInfo(object, type = c("I", "II", "Control", "I-Green", "I-Red"))
getManifestInfo(object, type = c("nLoci", "locusNames"))
getControlAddress(object, controlType = c("NORM_A", "NORM_C", "NORM_G", "NORM_T"))
}
\arguments{
  \item{object}{An object of either class \code{"RGChannelSet"} or
    \code{"RGChannelSetExtended"} or
    \code{"IlluminaMethylationManifest"}.  In case a data object is
    given (the two former possibilities), the associated manifest object
    is retrieved and used to obtain the information.}
  \item{type}{A single character describing what kind of information
    should be returned.  For \code{getProbeInfo} it represents the following
    subtypes of probes on the array: Type I, Type II, Controls as well
    as Type I (methylation measured in the Green channel) and Type II
    (methylation measured in the Red channel).  For
    \code{getManifestInfo} it represents either the number of
    methylation loci (approx. number of CpGs) on the array or the locus names.}
  \item{controlType}{A character vector of control types.}
}
\value{

  \code{getProbeData} returns the data slot of the manifest object
(mostly for internal use).  \code{getProbeInfo} returns a
\code{"data.frame"}, \code{getManifestInfo} returns either a single
number of a character vector and \code{getControlAddress} returns a
vector of addresses (probe locations).

}
\author{
Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\examples{
if (require(minfiData)) {

info <- getProbeInfo(RGsetEx, type = c("I"))
head(info)
info <- getProbeInfo(RGsetEx, type = c("Control"))
head(info)

}
}
