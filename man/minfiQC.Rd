\name{minfiQC}
\alias{minfiQC}
\title{
  easy one-step QC of methylation object
}
\description{
  This function combines a number of functions into a simple to use, one
  step QC step/
}
\usage{
minfiQC(object, fixOutliers = TRUE, verbose = FALSE)
}
\arguments{
  \item{object}{An object of class \code{[Genomic]MethylSet}.}
  \item{fixOutliers}{Should the function fix outlying observations
    (using \code{fixMethOutliers}) before running QC?}
  \item{verbose}{Should the function be verbose?}
}
\details{
  A number of functions are run sequentially on the \code{object}.

  First outlier values are thresholded using \code{fixMethOutliers}.
  Then qc is performed using \code{getQC} and then sample specific sex
  is estimated using \code{getSex}.
}
\value{
  A list with two values,
  \item{object}{The \code{object} processed by \code{fixMethOutliers}
    and with a column \code{predictedSex} added to the pheno data.}
  \item{qc}{A \code{DataFrame} with columns from the output of
  \code{getQC} and \code{getSex}}
}
\author{
  Kasper D. Hansen
}
\seealso{
  \code{\link{getSex}}, \code{\link{getQC}}, \code{\link{fixMethOutliers}}
}
\examples{
if(require(minfiData)) {
  out <- minfiQC(MsetEx)
  ## plotQC(out$qc)
  ## plotSex(out$sex)
}
}
