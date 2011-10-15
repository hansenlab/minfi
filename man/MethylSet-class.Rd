\name{MethylSet-class}
\Rdversion{1.1}
\docType{class}
\alias{MethylSet-class}
\alias{getMeth}
\alias{getUnmeth}
\alias{getBeta}
\alias{getM}
\alias{getAnno}
\alias{show,MethylSet-method}
\alias{initialize,MethylSet-method}
\title{Class MethylSet}
\description{
  This class holds preprocessed data for Illumina methylation microarrays.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("MethylSet", Meth, Unmeth, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{preprocessMethod}:}{Object of class \code{"character"}.
    This contains the preprocess method used to create the data.  This
    is a a character vector of length 3, the first component is the
    preprocessing method, followed by the version of \pkg{minfi} used to do
    the preprocessing followed by the version of the manifest package
    used to do the preprocessing.}
    \item{\code{assayData}:}{Object of class \code{"AssayData"}.  This
      has to contain the following two elements: \code{Meth} and \code{Unmeth}.}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}.}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}.}
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"}.}
    \item{\code{annotation}:}{Object of class \code{"character"}.
      Contains the type of array for this dataset.}
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"}.}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class \code{"eSet"}, distance 2.
Class \code{"\linkS4class{Versioned}"}, by class \code{"eSet"}, distance 3.
}
\section{Methods}{
  \describe{
    \item{\code{show}:}{The show method}
    \item{\code{initialize}:}{The initialize method, see above}
  } 
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\linkS4class{eSet}} for the basic class structure.  Objects of
  this class are for example created from an
  \code{\linkS4class{RGChannelSet}} using \code{\link{preprocessRaw}}.
}
\examples{
showClass("MethylSet")
}
\keyword{classes}
