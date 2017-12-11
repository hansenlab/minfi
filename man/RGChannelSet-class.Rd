\name{RGChannelSet-class}
\Rdversion{1.1}
\docType{class}
\alias{RGChannelSet-class}
\alias{RGChannelSetExtended-class}
\alias{RGChannelSet}
\alias{RGChannelSetExtended}
%% methods
\alias{show,RGChannelSet-method}
\alias{annotation,RGChannelSet-method}
\alias{annotation<-,RGChannelSet,ANY-method}
\alias{combine,RGChannelSet,RGChannelSet-method}
\alias{getBeta,RGChannelSet-method}
\alias{getManifest,RGChannelSet-method}
\alias{coerce,RGChannelSetExtended,RGChannelSet-method}
%% eSet-methods
\alias{pData,RGChannelSet-method}
\alias{pData<-,RGChannelSet,DataFrame-method}
\alias{sampleNames,RGChannelSet-method}
\alias{sampleNames<-,RGChannelSet,ANY-method}
\alias{featureNames,RGChannelSet-method}
\alias{featureNames<-,RGChannelSet-method}
%% functions
\alias{getGreen}
\alias{getRed}
\alias{getOOB}
\alias{getSnpBeta}
\alias{getNBeads}

\title{Class \code{"RGChannelSet"}}

\description{
  These classes represents raw (unprocessed) data from a two color
  micro array; specifically an Illumina methylation array.
}

\usage{
## Constructors

RGChannelSet(Green = new("matrix"), Red = new("matrix"),
             annotation = "", ...)

RGChannelSetExtended(Green = new("matrix"), Red = new("matrix"),
                    GreenSD = new("matrix"), RedSD = new("matrix"),
                    NBeads = new("matrix"), annotation = "", ...)

## Accessors

\S4method{annotation}{RGChannelSet}(object)
\S4method{pData}{RGChannelSet}(object)
\S4method{sampleNames}{RGChannelSet}(object)
\S4method{featureNames}{RGChannelSet}(object)
\S4method{getBeta}{RGChannelSet}(object, \dots)
getGreen(object)
getRed(object)
getNBeads(object)
\S4method{getManifest}{RGChannelSet}(object)

## Convenience functions
getOOB(object)
getSnpBeta(object)
}

\value{
  An object of class \code{RGChannelSet} or \code{RGChannelSetExtended}
  for the constructors.
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
  \item{annotation}{An annotation string, optional.}
  \item{\dots}{For the constructor(s), additional arguments to be passed
    to \code{SummarizedExperiment}; of particular interest are
    \code{colData}, \code{rowData} and \code{metadata}.  For
    \code{getBeta} these values gets passed onto \code{getBeta}.} 
}

\section{Constructors}{
  Instances are constructed using the
  \code{RGChannelSet} or \code{RGChannelSetExtended} functions with the
  arguments outlined above.

  \code{as(object, "RGChannelSet")} coerces a
  \code{RGChannelSetExtended} object into a \code{RGChannelSet}.
}

\section{Accessors}{
  \describe{
    \item{\code{getGreen}:}{Gets the Green channel as a matrix.}
    \item{\code{getRed}:}{Gets the Red channel as a matrix.}
    \item{\code{getNBeads}:}{Gets the number of beads as a matrix, this
      requires an \code{RGChannelSetExtended}.}
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
      around 3 distinct values corresponding to homo-, and
      hetero-zygotes.}
    \item{\code{combine}:}{Combines two different \code{RGChannelSet},
      eventually using the \code{combine} method for \code{eSet}.}
  }
}

\section{Tips}{
  The class inherits a number of useful methods from \code{SummarizedExperiment}.  In earliers versions of minfi, this class inherited
  from \code{eSet}, and we have kept of number of methods related to this, for example \code{pData}. 
  
  The best way to access phenotype data and sample names are \code{colData} and \code{colnames}.
  
  Amongst the useful methods are

  \describe{
    \item{\code{dim}, \code{nrow}, \code{ncol}}{The dimension (number of
    probes by number of samples) of the experiment.}
    \item{\code{colData}, \code{colnames}, \code{pData}, \code{sampleNames}}{Phenotype information and
    sample names.}
    \item{\code{rownames}, \code{featureNames}}{This is the addresses (probe identifiers)
    of the array.}
  }
}
  
\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}

\seealso{
  See \code{\linkS4class{SummarizedExperiment}} for the basic class that is used as a
  building block for \code{"RGChannelSet(Extended)"}.  See
  \code{\link{IlluminaMethylationManifest}} for a class representing the
  design of the array. 
}

\examples{
showClass("RGChannelSet")
}
