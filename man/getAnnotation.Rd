\name{getAnnotation}
\alias{getAnnotation}
\alias{getLocations}
\alias{getSnpInfo}
\alias{addSnpInfo}
\alias{getIslandStatus}
\alias{getProbeType}
\title{
  Accessing annotation for Illumina methylation objects
}
\description{
  These functions access provided annotation for various Illumina
  methylation objects.
}
\usage{
getAnnotation(object, what = "everything", lociNames = NULL,
              orderByLocation = FALSE, dropNonMapping = FALSE)

getLocations(object, mergeManifest = FALSE,
             orderByLocation = FALSE, lociNames = NULL)

getSnpInfo(object, snpType = NULL)
addSnpInfo(object, snpType = NULL)

getProbeType(object)

getIslandStatus(object, islandType = "UCSC")
}

\arguments{
  \item{object}{}
  \item{what}{}
  \item{lociNames}{}
  \item{orderByLocation}{}
  \item{dropNonMapping}{}
  \item{mergeManifest}{}
  \item{snpType}{}
  \item{islandType}{}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
For \code{getAnnotation}, a \code{DataFrame} with the requested
information.

For \code{getLocations}, a \code{GRanges} with the locations.

For \code{getProbeType} and \code{getIslandStatus}, a \code{character}
vector with the requested information.

For \code{getSnpInfo}, a \code{DataFrame} with the requested
information.  For \code{addSnpInfo}, an object of the same class as
\code{object} but with the SNP information added to the elementMetadata
of the granges of the object.
}

\author{
Kasper Daniel Hansen\email{khansen@jhsph.edu}
}
\seealso{
\code{\linkS4class{IlluminaMethylationAnnotation}}
}
\examples{
}