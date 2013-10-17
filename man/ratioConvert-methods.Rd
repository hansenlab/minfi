\name{ratioConvert-methods}
\docType{methods}
\alias{ratioConvert}
\alias{ratioConvert-methods}
\alias{ratioConvert,MethylSet-method}
\alias{ratioConvert,GenomicMethylSet-method}
\title{Converting methylation signals to ratios (Beta or M-values)}
\description{
  Converting methylation data from methylation and unmethylation
  channels, to ratios (Beta and M-values).
}

\usage{

\S4method{ratioConvert}{MethylSet}(object, what = c("beta", "M", "both"), keepCN = TRUE, \dots)
\S4method{ratioConvert}{GenomicMethylSet}(object, what = c("beta", "M", "both"), keepCN = TRUE, \dots)
}

\arguments{
  \item{object}{Either a \code{MethylSet}, or a \code{GenomicRatioSet}.}
  \item{what}{Which ratios should be computed and stored?}
  \item{keepCN}{A logical, should compy number values be computed and stored in the
    object?}
  \item{\dots}{Passed to \code{getBeta}, \code{getM} methods.}
}

\author{Kasper Daniel Hansen \email{khansen@jhsph.edu}}

\seealso{

  \code{\link{RatioSet}} or code{\link{GenomicRatioSet}} for the output
  object and \code{\link{MethylSet}} or code{\link{GenomicMethylSet}}
  for the input object.

}

\keyword{methods}

