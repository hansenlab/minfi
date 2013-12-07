\name{IlluminaMethylationManifest-class}
\Rdversion{1.1}
\docType{class}
\alias{IlluminaMethylationManifest-class}
\alias{IlluminaMethylationManifest}
\alias{show,IlluminaMethylationManifest-method}
\alias{getManifest,IlluminaMethylationManifest-method}
\alias{getManifest,character-method}
\alias{getManifest}
\alias{getProbeInfo}
\alias{getManifestInfo}
\alias{getControlAddress}
\alias{getControlTypes}

\title{Class \code{"IlluminaMethylationManifest"}}

\description{
  This is a class for representing an Illumina methylation microarray
  design, ie. the physical location and the probe sequences.  This
  information should be independent of genome build and annotation.
}

\usage{
## Constructor

IlluminaMethylationManifest(TypeI = new("DataFrame"),
                            TypeII = new("DataFrame"),
                            TypeControl = new("DataFrame"),
                            TypeSnpI = new("DataFrame"),
                            TypeSnpII = new("DataFrame"),
                            annotation = "")

## Data extraction
\S4method{getManifest}{IlluminaMethylationManifest}(object)
\S4method{getManifest}{character}(object)
getProbeInfo(object, type = c("I", "II", "Control",
                              "I-Green", "I-Red", "SnpI", "SnpII"))
getManifestInfo(object, type = c("nLoci", "locusNames"))
getControlAddress(object, controlType = c("NORM_A", "NORM_C",
                                          "NORM_G", "NORM_T"),
                  asList = FALSE)
getControlTypes(object)
}

\arguments{
  \item{object}{Either an object of class
    \code{IlluminaMethylationManifest} or class \code{character} for
  \code{getManifest}.  For \code{getProbeInfo}, \code{getManifestInfo}
  and \code{getControlAddress} an object of either class
  \code{RGChannelSet}, \code{IlluminaMethylationManifest}.}
  \item{TypeI}{A \code{DataFrame} of type I probes.}
  \item{TypeII}{A \code{DataFrame} of type II probes.}
  \item{TypeControl}{A \code{DataFrame} of control probes.}
  \item{TypeSnpI}{A \code{DataFrame} of SNP type I probes.}
  \item{TypeSnpII}{A \code{DataFrame} of SNP type II probes.}
  \item{annotation}{An annotation \code{character}.}
  \item{type}{A single character describing what kind of information
    should be returned.  For \code{getProbeInfo} it represents the following
    subtypes of probes on the array: Type I, Type II, Controls as well
    as Type I (methylation measured in the Green channel) and Type II
    (methylation measured in the Red channel).  For
    \code{getManifestInfo} it represents either the number of
    methylation loci (approx. number of CpGs) on the array or the locus names.}
  \item{controlType}{A character vector of control types.}
  \item{asList}{If \code{TRUE} the return object is a list with one
  component for each \code{controlType}.}
}

\section{Details}{
  The \code{data} slot contains the following objects: \code{TypeI},
  \code{TypeII} and \code{TypeControl} which are all of class
  \code{data.frame}, describing the array design.

  Methylation loci of type I are measured using two different probes, in
  either the red or the green channel.  The columns \code{AddressA},
  \code{AddresB} describes the physical location of the two probes on
  the array (with \code{ProbeSeqA}, \code{ProbeSeqB} giving the probe
  sequences), and the column \code{Color} describes which color channel
  is used.

  Methylation loci of type II are measured using a single probe, but
  with two different color channels.  The methylation signal is always
  measured in the green channel.
}

\section{Utilities}{
  In the following code, \code{object} is a \code{IlluminaMethylationManifest}.

  \describe{
    \item{\code{getManifest(object)}}{Get the manifest object.}
    \item{\code{getProbeInfo(object)}}{Returns a \code{data.frame}
    giving the type I, type II or control probes.  It is also possible
    to get the type I probes measured in the Green or Red channel.}
   \item{\code{getManifestInfo(object)}}{Get some information about the
    manifest object (the chip design).}
   \item{\code{getControlAddress(object)}}{Get the control addresses for
    control probes of a certain type.}
   \item{\code{getControlTypes(object)}}{Returns the types and the
    numbers of control probes of each type.}
  }
}

\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}

\seealso{
  \code{\link{IlluminaMethylationAnnotation}} for annotation
  information for the array (information depending on a specific genome
  build).
}
\examples{
if(require(IlluminaHumanMethylation450kmanifest)) {

show(IlluminaHumanMethylation450kmanifest)
head(getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "I"))
head(IlluminaHumanMethylation450kmanifest@data$TypeI)
head(IlluminaHumanMethylation450kmanifest@data$TypeII)
head(IlluminaHumanMethylation450kmanifest@data$TypeControl)

}
}
\keyword{classes}
