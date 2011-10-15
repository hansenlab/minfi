\name{IlluminaMethylationManifest-class}
\Rdversion{1.1}
\docType{class}
\alias{IlluminaMethylationManifest-class}
\alias{manifestNew}
\alias{show,IlluminaMethylationManifest-method}
\title{Class \code{"IlluminaMethylationManifest"}}
\description{
  This is a class for representing an Illumina methylation microarray
  design, ie. the physical location and the probe sequences.
}
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{new("IlluminaMethylationManifest", ...)}, but the preferred way is
to use \code{manifestNew()}.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{This environment holds the design objects.  See details.}
    \item{\code{annotation}:}{A character vector giving the annotation name.}
  }
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
\section{Methods}{
  \describe{
    \item{\code{show}:}{The show method.}
  }
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\link{getProbeData}}
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
