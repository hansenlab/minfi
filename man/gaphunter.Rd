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
gaphunter(object,threshold=0.05,keepoutliers=FALSE,outcutoff=0.01,verbose=TRUE)
}

\arguments{
  \item{object}{An object of class GenomicRatioSet, GenomicMethylSet, or matrix.
	If one of the first two, code{\link{getBeta}} is used to calculate beta values. If a 
	matrix, must be one of beta values.}
  \item{threshold}{The difference in consecutive, ordered beta values
    that defines the presence of a gap signal. Defaults to 5 percent.}
  \item{keepoutliers}{Should outlier-driven gap signals be kept in the
    results? Defaults to \code{FALSE}}
  \item{outcutoff}{Value used to identify gap signals driven by outliers.
    Defined as the percentage of the total sample size; the sum of
    samples in all groups except the largest must exceed this number of
    samples in order for the probe to still be considered a gap signal.
	Defaults to 1 percent.}
 \item{verbose}{logical value. If \code{TRUE}, it writes some messages
    indicating progress. If \code{FALSE} nothing should be printed.}
}
\details{
  
  The function can calculate a beta matrix or utilize a user-supplied matrix of
  beta values. 

  The function will idenfity probes with a gap in a beta signal greater than
  or equal to the defined \code{threshold}. This method provides an data-driven alternative
  to reference annotations in terms of accounting for polymorphism-affected probes.
  
  Outlier-driven gap signals are defined as those in which the sum of the smaller group(s) does
  not exceed a certain percentage of the sample size, defined by the argument outcutoff. 
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
  \item{algorithm}{A list detailing the arguments supplied to the function.}
}
\author{
Shan Andrews \email{sandre17@jhu.edu}.
}
\examples{
require(minfiData)
gmSet <- preprocessQuantile(MsetEx)
gapres <- gaphunter(gmSet,threshold=0.1,keepoutliers=TRUE)
}