\name{fixMethOutliers}
\alias{fixMethOutliers}
\title{
  Fix methylation outliers
}
\description{
  Methylation outliers (loci with very extreme values of the Meth or
  Unmeth channel) are identified and fixed (see details).
}
\usage{
fixMethOutliers(object, K = -3, verbose = FALSE)
}
\arguments{
  \item{object}{An object of class \code{[Genomic]MethylSet}.}
  \item{K}{The number of standard deviations away from the median when
    defining the outlier cutoff, see details.}
  \item{verbose}{Should the function be verbose?}
}
\details{
  This function fixes outlying methylation calls in the Meth channel and
  Unmeth channel separately.

  Unlike other types of arrays, all loci on a methylation array ought to
  measure something (apart from loci on the Y chromosome in a female
  sample).  An outlier is a loci with a very low value in one of the two
  methylation channels.  Typically, relatively few loci ought to be outliers.

  An outlier is defined in a sample and methylation channel specific
  way.  First the (sample, methylation channel) values are
  \code{log2(x+0.5)} transformed and then the median and mad of these
  values are computed.  An outlier is then defined to be any value less
  than the median plus \code{K} times the mad, and these outlier values
  are thresholded at the cutoff (on the original scale).
}
\value{
  An object of the same class as \code{object} where outlier values in
  the methylation channels have been thresholded.
}
\author{
  Rafael A. Irizarry and Kasper D. Hansen
}
\seealso{
  \code{\link{minfiQC}}
}
\examples{
if(require(minfiData)) {
  MsetEx <- fixMethOutliers(MsetEx)
}
}
