\name{getMethSignal}
\alias{getMethSignal}
\title{
  Various utilities
}
\description{
  Utility functions operating on objects from the minfi package.
}
\usage{
getMethSignal(object, what = c("Beta", "M"), \dots)
}
\arguments{
  \item{object}{An object from the minfi package supporting either
    \code{getBeta} or \code{getM}.}
  \item{what}{Which signal is returned.}
  \item{\dots}{Passed to the method described by argument \code{what}.}
}
\value{A matrix.}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\examples{
if(require(minfiData)) {
  getMethSignal(MsetEx, what = "Beta")
}
}
