\name{GenomicRatioSet-class}
\Rdversion{1.1}
\docType{class}
\alias{GenomicRatioSet-class}
\alias{GenomicRatioSet}
\alias{getCN}
%%methods
\alias{getBeta,GenomicRatioSet-method}
\alias{getM,GenomicRatioSet-method}
\alias{getCN,GenomicRatioSet-method}
\alias{show,GenomicRatioSet-method}
\alias{pData,GenomicRatioSet-method}
\alias{sampleNames,GenomicRatioSet-method}
\alias{featureNames,GenomicRatioSet-method}
\alias{annotation,GenomicRatioSet-method}
\alias{preprocessMethod,GenomicRatioSet-method}
\alias{mapToGenome,GenomicRatioSet-method}

\title{GenomicRatioSet instances}

\description{
  This class holds preprocessed data for Illumina methylation
  microarrays, mapped to a genomic location.
}

\usage{
## Constructor

GenomicRatioSet(gr = GRanges(), Beta = NULL, M = NULL,
                CN = NULL, pData = DataFrame(),
                annotation = "", preprocessMethod = "")

## Data extraction / Accessors

\S4method{getBeta}{GenomicRatioSet}(object)
\S4method{getM}{GenomicRatioSet}(object)
\S4method{getCN}{GenomicRatioSet}(object)
\S4method{pData}{GenomicRatioSet}(object)
\S4method{sampleNames}{GenomicRatioSet}(object)
\S4method{featureNames}{GenomicRatioSet}(object)
\S4method{annotation}{GenomicRatioSet}(object)
\S4method{preprocessMethod}{GenomicRatioSet}(object)
\S4method{mapToGenome}{GenomicRatioSet}(object, \dots)
}

\arguments{
  \item{object}{A \code{GenomicRatioSet}.}
  \item{gr}{A \code{GRanges} object.}
  \item{Beta}{A matrix of beta values (optional, see details).}
  \item{M}{A matrix of M values (optional, see details).}
  \item{CN}{A matrix of copy number values.}
  \item{pData}{A \code{DataFrame} or \code{data.frame} object.}
  \item{annotation}{An annotation character string.}
  \item{preprocessMethod}{A preprocess method character string.}
  \item{\dots}{For \code{mapToGenome}, this is ignored.}
}

\details{
  This class holds \code{M} or \code{Beta} values (or both) together
  with associated genomic coordinates.  It is not possible to get
  \code{Meth} or \code{Unmeth} values from this object.  The intention
  is to use this kind of object as an analysis end point.

  In case one of \code{M} or \code{Beta} is missing, the other is
  computed on the fly.  For example, M is computed from Beta as the
  logit (base 2) of the Beta values.
}

\section{Constructor}{
Instances are constructed using the \code{GenomicRatioSet} function with the
arguments outlined above.
}

\section{Accessors}{
  A number of useful accessors are inherited from
  \code{SummarizedExperiment}.
  
  In the following code, \code{object} is a \code{GenomicRatioSet}.
  
  \describe{
    \item{\code{getBeta(object)}}{Get Beta, see details.}
    \item{\code{getM(object)}}{get M-values, see details.}
    \item{\code{getCN(object)}}{get copy number, see details.}
    \item{\code{getManifest(object)}}{get the manifest associated with
      the object.}
    \item{\code{sampleNames(object)}, \code{featureNames(object)}}{Get
      the sampleNames (colnames) or the featureNames (rownames).}
    \item{\code{preprocessMethod(object)},
      \code{annotation(object)}}{Get the preprocess method or annotation
      \code{character}.}
  }
}

\section{Utilities}{
  \describe{
    \item{\code{mapToGenome(object)}}{ Since \code{object} is already
 	mapped to the genome, this method simply returns \code{object}
	unchanged.}
  }
}

\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}

\seealso{
  \code{\link{SummarizedExperiment}} for the basic class structure.
}

\examples{
  showClass("GenomicRatioSet")
}
