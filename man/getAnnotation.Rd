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

getSnpInfo(object, snpAnno = NULL)
addSnpInfo(object, snpAnno = NULL)

getProbeType(object)

getIslandStatus(object, islandAnno = NULL)
}

\arguments{
  \item{object}{A \code{minfi} object.}
  \item{what}{Which annotation objects should be returned?}
  \item{lociNames}{Restrict the return values to these loci.}
  \item{orderByLocation}{Should the return object be ordered according
    to genomic location.}
  \item{dropNonMapping}{Should loci that do not have a genomic location
    associated with it (by being marked as \code{unmapped} or
    \code{multi}) be dropped from the return object.}
  \item{mergeManifest}{Should the manifest be merged into the return
    object?}
  \item{snpAnno}{The snp annotation you want to use; \code{NULL}
    signifies picking the default.}
  \item{islandAnno}{Like \code{snpAnno}, but for islands.}
}

\details{

  \code{getAnnotation} returns requested annotation as a
  \code{DataFrame}, with each row corresponding to a methylation loci.
  If \code{object} is of class \code{IlluminaHumanAnnotation} no
  specific ordering of the return object is imposed.  If, on the other
  hand, the class of \code{object} imposes some natural order on the
  return object (ie. if the object is of class
  \code{[Genomic](Methyl|Ratio)Set}), this order is kept in the return
  object.  Note that \code{RGChannelSet} does not impose a specific
  ordering on the methylation loci.

  \code{getLocations} is a convenience function which returns
  \code{Locations} as a \code{GRanges} and which furthermore drops
  unmapped loci.  A user should not need to call this function, instead
  \code{mapToGenome} should be used to get genomic coordinates and
  \code{granges} to return those coordinates.

  To see which options are available for \code{what}, simply print the
  annotation object.
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
\code{\linkS4class{IlluminaMethylationAnnotation}} for the basic class,
\code{\link{mapToGenome}} for a better alternative (for users) to
\code{getLocations}. 
}
\examples{
if(require(minfiData)) {
  table(getIslandStatus(MsetEx))
  getAnnotation(MsetEx, what = "Manifest")
}
}