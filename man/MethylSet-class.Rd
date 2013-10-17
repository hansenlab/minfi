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
\alias{preprocessMethod}
%%methods
\alias{getMeth,MethylSet-method}
\alias{getUnmeth,MethylSet-method}
\alias{getBeta,MethylSet-method}
\alias{getM,MethylSet-method}
\alias{getCN,MethylSet-method}
\alias{preprocessMethod,MethylSet-method}
\alias{getManifest,MethylSet-method}
\alias{show,MethylSet-method}
\alias{initialize,MethylSet-method}
\alias{dropMethylationLoci}

\title{MethylSet instances}

\description{
  This class holds preprocessed data for Illumina methylation
  microarrays.
}

\usage{
## Constructor

MethylSet(Meth = new("matrix"), Unmeth = new("matrix"),  \dots)

## Data extraction / Accessors

\S4method{getMeth}{MethylSet}(object)
\S4method{getUnmeth}{MethylSet}(object)
\S4method{getBeta}{MethylSet}(object, type = "", offset = 0, betaThreshold = 0)
\S4method{getM}{MethylSet}(object, type = "", \dots)
\S4method{getCN}{MethylSet}(object, \dots)
\S4method{getManifest}{MethylSet}(object)
\S4method{preprocessMethod}{MethylSet}(object)

## Utilities
dropMethylationLoci(object, dropRS = TRUE, dropCH = TRUE)
}

\arguments{
  \item{object}{A \code{MethylSet}.}
  \item{Meth}{A matrix of methylation values (between zero and infinity)
    with each row being a methylation loci and each column a sample.}
  \item{Unmeth}{See the \code{Meth} argument.}
  \item{type}{How are the values calculated?  For \code{getBeta} setting
    \code{type="Illumina"} sets \code{offset=100} as per Genome Studio.
    For \code{getM} setting \code{type=""} computes M-values as the
    logarithm of \code{Meth}/\code{Unmeth}, otherwise it is computed as
    the logit of \code{getBeta(object)}.}
  \item{offset}{Offset in the beta ratio, see detail.}
  \item{betaThreshold}{Constrains the beta values to be in the inverval
    betwen \code{betaThreshold} and 1-\code{betaThreshold}.}  
  \item{dropRS}{Should SNP probes be dropped?}
  \item{dropCH}{Should CH probes be dropped}
  \item{\dots}{For the constructor: additional objects passes to the
    \code{eSet} constructor, particular a \code{phenoData} slot.  For
    \code{getM} these values gets passed onto \code{getBeta}.}
  
}

\details{
  This class inherits from \code{eSet}.  Essentially the class is a representation of a
  \code{Meth} matrix and a \code{Unmeth} matrix linked to a \code{pData} data frame.

  In addition, an annotation and a preprocessMethod slot is present.  The annotation slot describes
  the type of array and also which annotation package to use.  The preprocessMethod slot describes
  the kind of preprocessing that resulted in this dataset.

  A \code{MethylSet} stores \code{meth} and \code{Unmeth}.  From these it is easy to compute Beta
  values, defined as
  \deqn{\beta = \frac{\textrm{Meth}}{\textrm{Meth} + \textrm{Unmeth} + \textrm{offset}}}{%
    Beta = Meth / (Meth + Unmeth + offset)}
  The offset is chosen to avoid dividing with small values.  Illumina uses a default of 100.
  M-values (an unfortunate bad name) are defined as
  \deqn{M = \textrm{logit}(\beta) = \log(\textrm{Meth} / \textrm{Unmeth})}{%
    M = logit(Beta) = log( Meth / Unmeth )}
  This formula has problems if either Meth or Unmeth is zero.  For this reason, we can use
  \code{betaThreshold} to make sure Beta is neither 0 nor 1, before taken the logit.  What makes
  sense for the \code{offset} and \code{betaThreshold} depends crucially on how the data was
  preprocessed.  Do not expect the default values to be particular good.
}

\section{Constructor}{
Instances are constructed using the \code{MethylSet} function with the
arguments outlined above.
}

\section{Accessors}{
  In the following code, \code{object} is a \code{MethylSet}.

  \describe{
    \item{\code{getMeth(object)}, \code{getUnmeth(object)}}{Get the
      Meth or the Unmeth matrix }
    \item{\code{getBeta(object)}}{Get Beta, see details.}
    \item{\code{getM(object)}}{get M-values, see details.}
    \item{\code{getCN(object)}}{get copy number values which are defined
      as the sum of the methylation and unmethylation channel.}
    \item{\code{getManifest(object)}}{get the manifest associated with
      the object.}
    \item{\code{preprocessMethod(object)}}{Get the preprocess method \code{character}.}
      }
}

\section{Utilities}{
  In the following code, \code{object} is a \code{MethylSet}.
  
  \describe{
    \item{\code{dropMethylationLoci}(object)}{A unifed interface to removing methylation loci.  You
      can drop SNP probes (probes that measure SNPs, not probes containing SNPs) or CH probes (non-CpG
      methylation).}
  }
}

\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}

\seealso{
  \code{\linkS4class{eSet}} for the basic class structure.
  Objects of this class are typically created from an
  \code{\link{RGChannelSet}} using \code{\link{preprocessRaw}} or
  another preprocessing function.
}

\examples{
  showClass("MethylSet")
}

