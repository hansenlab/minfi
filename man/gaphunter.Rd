\name{gaphunter}
\alias{gaphunter}
\title{
  Find gap signals in 450k data
}
\description{
  This function finds probes in the Illumina 450k Array for which calculated beta values 
  cluster into distinct groups separated by a defined threshold. It identifies, for these
  \sQuote{gaps signals} the number of groups, the size of these groups, and the samples 
  in each group.
}
\usage{
gaphunter(object,Beta=NULL,threshold=0.05,keepoutliers=FALSE,outcutoff=0.01,verbose=TRUE)
}
\arguments{
  \item{object}{A \code{GenomicRatioSet}.}
  \item{Beta}{Matrix of beta values. If NULL, a beta matrix is calculated using code{\link{getBeta}}}
  \item{threshold}{The difference in consecutive, ordered beta values that defines the presence
	of a gap signal. Defaults to 5%.}
  \item{keepoutliers}{Should outlier-driven gap signal be kept in the results? Defaults to \code{FALSE}}
  \item{outcutoff}{Value used to identify gap signals driven by outliers. Defined as percentage
	of total samples size; the sum of individuals in all groups except the largest must exceed
	this number of people in order for the probe to still be considered a gap signal. Defaults
	to 1%}
  \item{verbose}{
    logical value. If \code{TRUE}, it writes out some messages
    indicating progress. If \code{FALSE} nothing should be printed.
  }
}
\details{
  
  There are two main ways to run the function. Either the \code{object} should
  be supplied and the \code{Beta} matrix is calculated, or the \code{Beta}
  can be supplied directly.

  The function will idenfity probes with a gap in a beta signal greater than
  or equal to the defined \code{threshold}. This method provides an data-driven alternative
  to reference annotations in terms of accounting for polymorphism-affected probes.
}
\value{
  A list with three values,
  \item{proberesults}{A data frame listing, for each identified gap signal, the number
	of groups and the size of each group.}
  \item{sampleresults}{a matrix of dimemsions probes (rows) by samples (columns). Individuals
	are assigned numbers based onthe groups into which they cluster. Lower number groups
	indicate lower mean methylation values for the group. For example, individuals coded as
	\sQuote{1} will have a lower mean methylation value than those individuals coded as 
	\sQuote{2}.}
}
\author{
Shan Andrews \email{sandre17@jhu.edu}.
}
\examples{
if(require(minfiData)) {
  gmSet <- preprocessQuantile(MsetEx)
  gapres <- gaphunter(gmSet)
}
}
