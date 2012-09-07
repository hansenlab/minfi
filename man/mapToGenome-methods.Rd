\name{mapToGenome-methods}
\docType{methods}
\alias{mapToGenome}
\alias{mapToGenome-methods}
\alias{mapToGenome,MethylSet-method}
\alias{mapToGenome,RGChannelSet-method}
\title{Mapping methylation data to the genome}
\description{
  Mapping Ilumina methylation array data to the genome using an
  annotation package.  Depending on the genome, not all methylation loci
  may have a genomic position.
}

\usage{

\S4method{mapToGenome}{MethylSet}(object, genomeBuild = c("hg19", "hg18"), drop = TRUE, mergeManifest = FALSE)
\S4method{mapToGenome}{RGChannelSet}(object, ...)
}

\arguments{
  \item{object}{Either a \code{MethylSet} or \code{RGChannelSet}.}
  \item{genomeBuild}{Which version of the genome to use?}
  \item{drop}{Should unmapped loci be dropped?}
  \item{mergeManifest}{Should the information in the associated manifest
    package be merged into the location \code{GRanges}?}
  \item{...}{Passed to the method for \code{MethylSet}.}
}

\details{
  FIXME: details on the MethylSet method.
  
  The \code{RGChannelSet} method of this function is a convenience
  function: the \code{RGChannelSet} is first transformed into a
  \code{MethylSet} using \code{preprocessRaw}.  The resulting
  \code{MethylSet} is then mapped directly to the genome.
}

\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}

\seealso{
  \code{\link{GenomicMethylSet}} for the output object and
  \code{\link{MethylSet}} for the input object.  Also,
  \code{\link{getLocations}} obtains the genomic locations for a given object. 
}

\keyword{methods}

