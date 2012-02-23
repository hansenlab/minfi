\name{RGChannelSet-class}
\Rdversion{1.1}
\docType{class}
\alias{RGChannelSet-class}
\alias{RGChannelSetExtended-class}
\alias{RGChannelSet}
\alias{RGChannelSetExtended}
\alias{show,RGChannelSet-method}
\alias{initialize,RGChannelSet-method}
\alias{initialize,RGChannelSetExtended-method}
\alias{getGreen}
\alias{getRed}
\alias{getManifest}
\title{Class \code{"RGChannelSet"}}
\description{
  These classes represents raw (unprocessed) data from a two color
  micro array; specifically an Illumina methylation array.
}
\section{Objects from the Class}{
Objects can be created by the \code{RGChannelSet} or
\code{RGChannelSetExtended} constructors, ie. by calls of the form
\code{RGChannelSet(Red, Green, ...)} or \code{RGChannelSetExtended(Red,
  Green, RedSD, GreenSD, NBeads, ...)}. 
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"}.  A
      \code{RGChannelSet} contains elements with names \code{Green} and
      \code{Red} and a \code{"RGChannelSetExtended"} additionally adds
      \code{GreenSD}, \code{RedSD} and \code{NBeads}.}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}.}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}.}
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"}.}
    \item{\code{annotation}:}{Object of class \code{"character"}.
      Contains the array name (without \code{"manifest"}.)}
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"}.}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class \code{"eSet"}, distance 2.
Class \code{"\linkS4class{Versioned}"}, by class \code{"eSet"}, distance 3.
}
\usage{
RGChannelSet(Green = new("matrix"), Red = new("matrix"), ...)
RGChannelSetExtended(Green = new("matrix"), Red = new("matrix"),
                     GreenSD = new("matrix"), RedSD = new("matrix"),
                     NBeads = new("matrix"), ...)
}
\arguments{
  \item{Green}{A matrix of Green channel values (between zero and infinity)
    with each row being a methylation loci and each column a sample.}
  \item{Red}{See the \code{Green} argument, but for the Green channel.}
  \item{GreenSD}{See the \code{Green} argument, but for standard
    deviations of the Green channel summaries.}  
  \item{RedSD}{See the \code{Green}, but for standard deviations of the
    Red channel summaries.} 
  \item{NBeads}{See the \code{Green} argument, but contains the number
    of beads used to summarize the Green and Red channels.}
  \item{...}{Additional objects passes to the \code{eSet} constructor,
    particular a \code{phenoData} slot.}
}
\section{Methods}{
  \describe{
    \item{\code{show}:}{The show method}
    \item{\code{initialize}:}{The initialize method, see above}
  }
}
\section{Utilities}{
  \describe{
    \item{\code{getGreen}:}{Gets the Green channel as a matrix.}
    \item{\code{getRed}:}{Gets the Red channel as a matrix.}
    \item{\code{getManifest}:}{Gets the manifest object associated with
      the array type}
  }
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  See \code{\linkS4class{eSet}} for the basic class that is used as a
  building block for \code{"RGChannelSet(Extended)"}.  See
  \code{\linkS4class{IlluminaMethylationManifest}} for a class
  representing the design of the array.
}
\examples{
showClass("RGChannelSet")
}
\keyword{classes}
