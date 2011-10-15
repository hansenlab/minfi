\name{detectionP}
\alias{detectionP}
\title{
  Detection p-values for all probed genomic positions.
}
\description{
  This function identifies failed positions defined as both the methylated and unmethylated channel reporting background signal levels. 
}
\usage{
detectionP(rgSet, type = "m+u")
}
\arguments{
  \item{rgSet}{An \code{RGChannelSet}.}
  \item{type}{How to calculate p-values. Only \code{m+u} is currently implemented (See details).}
}
\details{
  
  A detection p-value is returned for every genomic position in every sample. Small
  p-values indicate a good position. Positions with non-significant p-values
  (typically >0.01) should not be trusted.

  The \code{m+u} method compares the total DNA signal (Methylated +
  Unmethylated) for each position to the background signal level. The
  background is estimated using negative control positions, assuming a
  normal distribution. Calculations are performed on the original
  (non-log) scale.

  This function is different from the detection routine in Genome Studio.
}
\value{
A matrix with detection p-values.
}
\author{
Martin Aryee \email{aryee@jhu.edu}.
}
\examples{
if (require(minfiData)) {

detP <- detectionP(RGsetEx)
failed <- detP>0.01
colMeans(failed) # Fraction of failed positions per sample
sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?

}
}
