\name{MethylSet-class}
\Rdversion{1.1}
\docType{class}
\alias{MethylSet-class}
\alias{MethylSet}
%% generics
\alias{getMeth}
\alias{getUnmeth}
\alias{getBeta}
\alias{getM}
%%methods
\alias{getMeth,MethylSet-method}
\alias{getUnmeth,MethylSet-method}
\alias{getBeta,MethylSet-method}
\alias{getM,MethylSet-method}
\alias{show,MethylSet-method}
\alias{initialize,MethylSet-method}

\title{MethylSet instances}

\description{
  This class holds preprocessed data for Illumina methylation
  microarrays.
}

\usage{
## Constructor
MethylSet(Meth = new("matrix"), Unmeth = new("matrix"),  ...)

## Data extraction / Accessors
\S4method{getMeth}{MethylSet}(object)
\S4method{getUnmeth}{MethylSet}(object)
\S4method{getBeta}{MethylSet}(object, type = c("minfi", "Illumina"))
\S4method{getM}{MethylSet}(object, type = c("minfi", "Illumina"))
}

\arguments{
  \item{object}{A \code{MethylSet}.}
  \item{Meth}{A matrix of methylation values (between zero and infinity)
    with each row being a methylation loci and each column a sample.}
  \item{Unmeth}{See the \code{Meth} argument.}
  \item{...}{Additional objects passes to the \code{eSet} constructor,
    particular a \code{phenoData} slot.}
  \item{type}{How is the values calculated?}
}

\details{
  FIXME: Mention the preprocess, annotation slots
}

\section{Constructor}{
Instances are constructed using the \code{MethylSet} function with the
arguments outlined above.
}

\section{Accessors}{
  In the following code, \code{object} is a \code{MethylSet}.

  \describe{
    \item{\code{getMeth(object)}, \code{getUnmeth(object)}}{Get the
    }
    \item{\code{getBeta(object)}}{}
    \item{\code{getM(object)}}{}
  }
}

\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}

\seealso{
  \code{\linkS4class{eSet}} for the basic class structure.
  Objects of this class are for example created from an
  \code{\link{RGChannelSet}} using \code{\link{preprocessRaw}}.
}

\examples{
  showClass("MethylSet")
}
