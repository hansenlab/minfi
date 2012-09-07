\name{GenomicMethylSet-class}
\Rdversion{1.1}
\docType{class}
\alias{GenomicMethylSet-class}
\alias{GenomicMethylSet}
%%methods
\alias{getMeth,GenomicMethylSet-method}
\alias{getUnmeth,GenomicMethylSet-method}
\alias{getBeta,GenomicMethylSet-method}
\alias{getM,GenomicMethylSet-method}
\alias{show,GenomicMethylSet-method}
\alias{pData,GenomicMethylSet-method}
\alias{sampleNames,GenomicMethylSet-method}
\alias{featureNames,GenomicMethylSet-method}
\alias{annotation,GenomicMethylSet-method}
\alias{preprocessMethod,GenomicMethylSet-method}
\alias{granges,GenomicMethylSet-method}
\alias{start,GenomicMethylSet-method}
\alias{end,GenomicMethylSet-method}
\alias{width,GenomicMethylSet-method}
\alias{strand,GenomicMethylSet-method}
\alias{seqnames,GenomicMethylSet-method}
\alias{seqlevels,GenomicMethylSet-method}
\alias{seqlengths,GenomicMethylSet-method}
\alias{genome,GenomicMethylSet-method}

\title{GenomicMethylSet instances}

\description{
  This class holds preprocessed data for Illumina methylation
  microarrays, mapped to a genomic location.
}

\usage{
## Constructor

GenomicMethylSet(gr, Meth, Unmeth, pData, annotation, preprocessMethod)

## Data extraction / Accessors

\S4method{getMeth}{GenomicMethylSet}(object)
\S4method{getUnmeth}{GenomicMethylSet}(object)
\S4method{getBeta}{GenomicMethylSet}(object, type = "", offset = 0, betaThreshold = 0)
\S4method{getM}{GenomicMethylSet}(object, ...)
\S4method{pData}{GenomicMethylSet}(object)
\S4method{sampleNames}{GenomicMethylSet}(object)
\S4method{featureNames}{GenomicMethylSet}(object)
\S4method{annotation}{GenomicMethylSet}(object)
\S4method{preprocessMethod}{GenomicMethylSet}(object)
\S4method{granges}{GenomicMethylSet}(x, ...)
\S4method{start}{GenomicMethylSet}(x, ...)
\S4method{end}{GenomicMethylSet}(x, ...)
\S4method{width}{GenomicMethylSet}(x)
\S4method{strand}{GenomicMethylSet}(x, ...)
\S4method{seqnames}{GenomicMethylSet}(x)
\S4method{seqlevels}{GenomicMethylSet}(x)
\S4method{seqlengths}{GenomicMethylSet}(x)
\S4method{genome}{GenomicMethylSet}(x)
}

\arguments{
  \item{object}{A \code{GenomicMethylSet}.}
  \item{x}{A \code{GenomicMethylSet}.}
  \item{gr}{A \code{GRanges} object.}
  \item{pData}{A \code{DataFrame} or \code{data.frame} object.}
  \item{Meth}{A matrix of methylation values (between zero and infinity)
    with each row being a methylation loci and each column a sample.}
  \item{Unmeth}{See the \code{Meth} argument.}
  \item{annotation}{An annotation character string.}
  \item{preprocessMethod}{A preprocess method character string.}
  \item{type}{How are the values calculated?  For \code{getBeta} setting
  \code{type="Illumina"} sets \code{offset=100} as per Genome Studio.
  For \code{getM} setting \code{type=""} computes M-values as the
  logarithm of \code{Meth}/\code{Unmeth}, otherwise it is computed as
  the logit of \code{getBeta(object)}.}
  \item{offset}{Offset in the beta ratio, see detail.}
  \item{betaThreshold}{Constrains the beta values to be in the inverval
    betwen \code{betaThreshold} and 1-\code{betaThreshold}.}  
  \item{...}{ For \code{getM} these values gets passed onto
  \code{getBeta}.} 
}

\details{
  For a detailed discussion of \code{getBeta} and \code{getM} see the
  deatils section of \code{\link{MethylSet}}.
}

\section{Constructor}{
Instances are constructed using the \code{GenomicMethylSet} function with the
arguments outlined above.
}

\section{Accessors}{
  In the following code, \code{object} is a \code{GenomicMethylSet}.
  
  \describe{
    \item{\code{getMeth(object)}, \code{getUnmeth(object)}}{Get the
      Meth or Unmeth matrix.}
    \item{\code{getBeta(object)}}{Get Beta, see details.}
    \item{\code{getM(object)}}{get M-values, see details.}
    \item{\code{getManifest(object)}}{get the manifest associated with
      the object.}
    \item{\code{sampleNames(object)}, \code{featureNames(object)}}{Get
      the sampleNames (colnames) or the featureNames (rownames).}
    \item{\code{preprocessMethod(object)},
      \code{annotation(object)}}{Get the preprocess method or annotation
      \code{character}.}
    \item{\code{granges(object)}}{Get the genomic locations associated
      with the object.}
    \item{\code{start}, \code{end}, \code{width}, \code{strand},
      \code{seqnames}, \code{seqlevels}, \code{seqlengths},
    \code{genome}}{Works directly on \code{granges(object)}.}
  }
}

\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}

\seealso{
  \code{\link{SummarizedExperiment}} for the basic class structure.
  Objects of this class are typically created by using the function
  \code{\link{mapToGenome}} on a \code{\link{MethylSet}}.}

\examples{
  showClass("GenomicMethylSet")
}
