\name{RGChannelSet-class}
\Rdversion{1.1}
\docType{class}
\alias{RGChannelSet-class}
\alias{RGChannelSetExtended-class}
\alias{RGChannelSet}
\alias{RGChannelSetExtended}
%% methods
\alias{show,RGChannelSet-method}
\alias{initialize,RGChannelSet-method}
\alias{initialize,RGChannelSetExtended-method}
\alias{getBeta,RGChannelSet-method}
\alias{getManifest,RGChannelSet-method}
%% functions
\alias{getGreen}
\alias{getRed}
\alias{getOOB}
\alias{getSnpBeta}

\title{Class \code{"RGChannelSet"}}

\description{
  These classes represents raw (unprocessed) data from a two color
  micro array; specifically an Illumina methylation array.
}

\usage{
## Constructors

RGChannelSet(Green = new("matrix"), Red = new("matrix"), \dots)
RGChannelSetExtended(Green = new("matrix"), Red = new("matrix"),
                     GreenSD = new("matrix"), RedSD = new("matrix"),
                     NBeads = new("matrix"), \dots)

## Accessors

\S4method{getBeta}{RGChannelSet}(object, \dots)
getGreen(object)
getRed(object)
\S4method{getManifest}{RGChannelSet}(object)

## Convenience functions
getOOB(object)
getSnpBeta(object)
}

\arguments{
  \item{object}{An \code{RGChannelSet} (or \code{RGChannelSetExtended}).}
  \item{Green}{A matrix of Green channel values (between zero and infinity)
    with each row being a methylation loci and each column a sample.}
  \item{Red}{See the \code{Green} argument, but for the Green channel.}
  \item{GreenSD}{See the \code{Green} argument, but for standard
    deviations of the Green channel summaries.}  
  \item{RedSD}{See the \code{Green}, but for standard deviations of the
    Red channel summaries.} 
  \item{NBeads}{See the \code{Green} argument, but contains the number
    of beads used to summarize the Green and Red channels.}
  \item{\dots}{Additional objects passes to the \code{eSet} constructor,
    particular a \code{phenoData} slot.}
}

\section{Constructors}{
  Instances are constructed using the
  \code{RGChannelSet} or \code{RGChannelSetExtended} functions with the
  arguments outlined above.
}

\section{Accessors}{
  \describe{
    \item{\code{getGreen}:}{Gets the Green channel as a matrix.}
    \item{\code{getRed}:}{Gets the Red channel as a matrix.}
    \item{\code{getManifest}:}{Gets the manifest object itself
      associated with the array type}
  }
}

\section{Convenience functions}{
  \describe{
    \item{\code{getOOB}:}{Retrives the so-called \dQuote{out-of-band}
      (OOB) probes.  These are the measurements of Type I probes in the
      \dQuote{wrong} color channel.  Return value is a list with two
      matrices, named \code{Red} and \code{Grn}.}
    \item{\code{getSnpBeta}:}{Retrives the measurements of the 65 SNP
      probes located on the array.  These SNP probes are intended to be used
      for sample tracking and sample mixups.  The return value is a
      matrix of beta values.  Each SNP probe ought to have values clustered
      around 3 distinct values corresponding to homo-, and hetero-zygotes.}
  }
}

\section{Tips}{
  The class inherits a number of useful methods from \code{eSet}.
  Amongst these are

  \describe{
    \item{\code{dim}, \code{nrow}, \code{ncol}}{The dimension (number of
    probes by number of samples) of the experiment.}
    \item{\code{pData}, \code{sampleNames}}{Phenotype information and
    sample names.}
    \item{\code{featureNames}}{This is the addresses (probe identifiers)
    of the array.}
  }
}
  
\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}

\seealso{
  See \code{\linkS4class{eSet}} for the basic class that is used as a
  building block for \code{"RGChannelSet(Extended)"}.  See
  \code{\link{IlluminaMethylationManifest}} for a class representing the
  design of the array. 
}

\examples{
showClass("RGChannelSet")
}
