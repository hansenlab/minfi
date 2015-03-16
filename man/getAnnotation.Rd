\name{getAnnotation}
\alias{getAnnotation}
\alias{getAnnotationObject}
\alias{getLocations}
\alias{getSnpInfo}
\alias{addSnpInfo}
\alias{dropLociWithSnps}
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

getAnnotationObject(object)

getSnpInfo(object, snpAnno = NULL)
addSnpInfo(object, snpAnno = NULL)

dropLociWithSnps(object, snps = c("CpG", "SBE"), maf = 0, snpAnno = NULL)

getProbeType(object, withColor = FALSE)

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
  \item{withColor}{Should the return object have the type I probe color
    labelled?}
  \item{snps}{The type of SNPs used.}
  \item{maf}{Minor allelle fraction.}
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

  \code{getAnnotationObject} returns the annotation object, as opposed
  to the annotation the object contains.  This is useful for printing
  and examining the contents of the object.

  \code{getLocations} is a convenience function which returns
  \code{Locations} as a \code{GRanges} and which furthermore drops
  unmapped loci.  A user should not need to call this function, instead
  \code{mapToGenome} should be used to get genomic coordinates and
  \code{granges} to return those coordinates.

  \code{getSnpInfo} is a conevnience function which gets a SNP DataFrame
  containing information on which probes contains SNPs
  where. \code{addSnpInfo} adds this information to the \code{rowRanges}
  or \code{granges} of the object.  \code{dropLociWithSnps} is a
  convenience function for removing loci with SNPs based on their MAF.
  
  To see which options are available for \code{what}, simply print the
  annotation object, possibly using \code{getAnnotationObject}.
}

\value{

  For \code{getAnnotation}, a \code{DataFrame} with the requested
  information.

  For \code{getAnnotationObject}, a \code{IlluminaMethylationAnnotation}
  object. 

  For \code{getLocations}, a \code{GRanges} with the locations.

  For \code{getProbeType} and \code{getIslandStatus}, a \code{character}
  vector with the requested information.
  
  For \code{getSnpInfo}, a \code{DataFrame} with the requested
  information.  For \code{addSnpInfo}, an object of the same class as
  \code{object} but with the SNP information added to the elementMetadata
  of the granges of the object.

  For \code{dropLociWithSnps} an object of the same kind as the input,
  possibly with fewer loci.
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
