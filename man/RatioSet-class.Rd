\name{RatioSet-class}
\Rdversion{1.1}
\docType{class}
\alias{RatioSet-class}
\alias{RatioSet}
%%methods
\alias{show,RatioSet-method}
\alias{preprocessMethod,RatioSet-method}
\alias{annotation,RatioSet-method}
\alias{annotation<-,RatioSet,ANY-method}
\alias{getBeta,RatioSet-method}
\alias{getM,RatioSet-method}
\alias{getCN,RatioSet-method}
\alias{combine,RatioSet,RatioSet-method}
%% eSet-methods
\alias{pData,RatioSet-method}
\alias{pData<-,RatioSet,DataFrame-method}
\alias{sampleNames,RatioSet-method}
\alias{sampleNames<-,RatioSet,ANY-method}
\alias{featureNames,RatioSet-method}
\alias{featureNames<-,RatioSet-method}

\title{RatioSet instances}

\description{
  This class holds preprocessed data for Illumina methylation
  microarrays.
}

\usage{
## Constructor

RatioSet(Beta = NULL, M = NULL, CN = NULL,
        annotation = "", preprocessMethod = "", ...)

## Data extraction / Accessors

\S4method{getBeta}{RatioSet}(object)
\S4method{getM}{RatioSet}(object)
\S4method{getCN}{RatioSet}(object)
\S4method{preprocessMethod}{RatioSet}(object)
\S4method{annotation}{RatioSet}(object)
\S4method{pData}{RatioSet}(object)
\S4method{sampleNames}{RatioSet}(object)
\S4method{featureNames}{RatioSet}(object)
}

\arguments{
  \item{object}{A \code{RatioSet}.}
  \item{Beta}{A matrix of beta values (between zero and one)
    with each row being a methylation loci and each column a sample.}
  \item{M}{A matrix of log-ratios (between minus infinity and infinity)
    with each row being a methylation loci and each column a sample.}
  \item{CN}{An optional matrix of copy number estimates
    with each row being a methylation loci and each column a sample.}
    \item{annotation}{An annotation string, optional.}
  \item{preprocessMethod}{A \code{character}, optional.}
  \item{\dots}{For the constructor, additional arguments to be passed to
    \code{SummarizedExperiment}; of particular interest are
    \code{colData}, \code{rowData} and \code{metadata}.
    For \code{getM} these values gets passed onto
    \code{getBeta}.} 
}

\value{
  An object of class \code{RatioSet} for the constructor.
}

\details{
  This class inherits from \code{eSet}.  Essentially the class
  is a representation of a \code{Beta} matrix and/or a \code{M} matrix
  and optionally a \code{CN} (copy number) matrix linked to a
  \code{pData} data frame.

  In addition, an annotation and a preprocessMethod slot is present.  The annotation slot describes
  the type of array and also which annotation package to use.  The preprocessMethod slot describes
  the kind of preprocessing that resulted in this dataset.

  For a \code{RatioSet}, M-values are defined as \code{logit2} of the
  Beta-values if the M-values are not present in the object.  Similarly,
  if only M-values are present in the object, Beta-values are
  \code{ilogit2} of the M-values.
}

\section{Constructor}{
Instances are constructed using the \code{RatioSet} function with the
arguments outlined above.
}

\section{Accessors}{
  In the following code, \code{object} is a \code{RatioSet}.

  \describe{
    \item{\code{getBeta(object)}, \code{getM(object)}, \code{CN(object)}}{Get the
      Beta, M or CN matrix.}
    \item{\code{getManifest(object)}}{get the manifest associated with
      the object.}
    \item{\code{preprocessMethod(object)}}{Get the preprocess method \code{character}.}
      }
}


\section{Utilities}{
  In the following code, \code{object} is a \code{RatioSet}.
  
  \describe{
    \item{\code{combine}:}{Combines two different \code{RatioSet},
      eventually using the \code{combine} method for \code{eSet}.}
  }
}

\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}

\seealso{
  \code{\linkS4class{eSet}} for the basic class structure.
  Objects of this class are typically created from an
  \code{\link{MethylSet}} using \code{\link{ratioConvert}}.
}

\examples{
  showClass("RatioSet")
}

