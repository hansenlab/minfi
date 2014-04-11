\name{mapToGenome-methods}
\docType{methods}
\alias{mapToGenome}
\alias{mapToGenome-methods}
\alias{mapToGenome,MethylSet-method}
\alias{mapToGenome,RatioSet-method}
\alias{mapToGenome,RGChannelSet-method}
\title{Mapping methylation data to the genome}
\description{
  Mapping Ilumina methylation array data to the genome using an
  annotation package.  Depending on the genome, not all methylation loci
  may have a genomic position.
}

\usage{

\S4method{mapToGenome}{MethylSet}(object, mergeManifest = FALSE)
\S4method{mapToGenome}{MethylSet}(object, mergeManifest = FALSE)
\S4method{mapToGenome}{RGChannelSet}(object, \dots)
}

\arguments{
  \item{object}{Either a \code{MethylSet}, a \code{RGChannelSet} or a
    \code{RatioSet}.} 
  \item{mergeManifest}{Should the information in the associated manifest
    package be merged into the location \code{GRanges}?}
  \item{\dots}{Passed to the method for \code{MethylSet}.}
}

\details{
  FIXME: details on the MethylSet method.
  
  The \code{RGChannelSet} method of this function is a convenience
  function: the \code{RGChannelSet} is first transformed into a
  \code{MethylSet} using \code{preprocessRaw}.  The resulting
  \code{MethylSet} is then mapped directly to the genome.

  This function silently drops loci which cannot be mapped to a genomic
  position, based on the associated annotation package.
}

\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}

\seealso{
  \code{\link{GenomicMethylSet}} for the output object and
  \code{\link{MethylSet}} for the input object.  Also,
  \code{\link{getLocations}} obtains the genomic locations for a given object. 
}

\keyword{methods}

