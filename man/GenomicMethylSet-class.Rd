\name{GenomicMethylSet-class}
\Rdversion{1.1}
\docType{class}
\alias{GenomicMethylSet-class}
\alias{GenomicMethylSet}
%%methods
\alias{show,GenomicMethylSet-method}
\alias{annotation,GenomicMethylSet-method}
\alias{annotation<-,GenomicMethylSet,ANY-method}
\alias{preprocessMethod,GenomicMethylSet-method}
\alias{getMeth,GenomicMethylSet-method}
\alias{getUnmeth,GenomicMethylSet-method}
\alias{getBeta,GenomicMethylSet-method}
\alias{getM,GenomicMethylSet-method}
\alias{getCN,GenomicMethylSet-method}
\alias{mapToGenome,GenomicMethylSet-method}
\alias{combine,GenomicMethylSet,GenomicMethylSet-method}
%% eSet-methods
\alias{pData,GenomicMethylSet-method}
\alias{pData<-,GenomicMethylSet,DataFrame-method}
\alias{sampleNames,GenomicMethylSet-method}
\alias{sampleNames<-,GenomicMethylSet,ANY-method}
\alias{featureNames,GenomicMethylSet-method}
\alias{featureNames<-,GenomicMethylSet-method}

\title{GenomicMethylSet instances}

\description{
  This class holds preprocessed data for Illumina methylation
  microarrays, mapped to a genomic location.
}

\usage{
## Constructor

GenomicMethylSet(gr = GRanges(), Meth = new("matrix"),
                 Unmeth = new("matrix"), annotation = "",
                 preprocessMethod = "", ...)

## Data extraction / Accessors

\S4method{getMeth}{GenomicMethylSet}(object)
\S4method{getUnmeth}{GenomicMethylSet}(object)
\S4method{getBeta}{GenomicMethylSet}(object, type = "", offset = 0, betaThreshold = 0)
\S4method{getM}{GenomicMethylSet}(object, type = "", \dots)
\S4method{getCN}{GenomicMethylSet}(object, \dots)
\S4method{pData}{GenomicMethylSet}(object)
\S4method{sampleNames}{GenomicMethylSet}(object)
\S4method{featureNames}{GenomicMethylSet}(object)
\S4method{annotation}{GenomicMethylSet}(object)
\S4method{preprocessMethod}{GenomicMethylSet}(object)
\S4method{mapToGenome}{GenomicMethylSet}(object, \dots)
}

\arguments{
  \item{object}{A \code{GenomicMethylSet}.}
  \item{gr}{A \code{GRanges} object.}
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
  \item{\dots}{For the constructor, additional arguments to be passed to
    \code{SummarizedExperiment}; of particular interest are
    \code{colData} and \code{metadata}. For \code{getM} these values
    gets passed onto \code{getBeta}. For \code{mapToGenome}, this is
    ignored.}
}

\value{
  An object of class \code{GenomicMethylSet} for the constructor.
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
  A number of useful accessors are inherited from
  \code{RangedSummarizedExperiment}.
  
  In the following code, \code{object} is a \code{GenomicMethylSet}.
  
  \describe{
    \item{\code{getMeth(object)}, \code{getUnmeth(object)}}{Get the
      Meth or Unmeth matrix.}
    \item{\code{getBeta(object)}}{Get Beta, see details.}
    \item{\code{getM(object)}}{get M-values, see details.}
    \item{\code{getCN(object)}}{get copy number values which are defined
      as the sum of the methylation and unmethylation channel.}
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
    \item{\code{combine}:}{Combines two different \code{GenomicMethylSet},
      eventually using the \code{cbind} method for \code{SummarizedExperiment}.}
  }
}


\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}

\seealso{
  \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} in the
  \pkg{SummarizedExperiment} package for the basic class structure.
  Objects of this class are typically created by using the function
  \code{\link{mapToGenome}} on a \code{\link{MethylSet}}.}

\examples{
  showClass("GenomicMethylSet")
}
