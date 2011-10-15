\name{preprocessIllumina}
\alias{preprocessIllumina}
\alias{bgcorrect.illumina}
\alias{normalize.illumina.control}
\title{
  Perform preprocessing as Genome Studio.
}
\description{
  These functions implements preprocessing for Illumina methylation
  microarrays as used in Genome Studio, the standard software provided
  by Illumina.
}
\usage{
preprocessIllumina(rgSet, bg.correct = TRUE, normalize = c("controls", "no"),
    reference = 1)
bgcorrect.illumina(rgSet)
normalize.illumina.control(rgSet, reference = 1)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet}.}
  \item{bg.correct}{logical, should background correction be performed?}
  \item{normalize}{logical, should (control) normalization be performed?}
  \item{reference}{for control normalization, which array is the reference?}
}
\details{
  We have reverse engineered the preprocessing methods from Genome
  Studio, based on the documentation.

  The current implementation of control normalization is equal to what
  Genome Studio provides (this statement is based on comparing Genome
  Studio output to the output of this function), with the following
  caveat: this kind of normalization requires the selection of a
  reference array.  It is unclear how Genome Studio selects the
  reference array, but we allow for the manual specification of this
  parameter.

  The current implementation of background correction is roughly equal
  to Genome Studio.  Based on examining the output of 24 arrays, we are
  able to exactly recreate 18 out of the 24.  The remaining 6 arrays had
  a max discrepancy in the Red and/or Green channel of 1-4 (this is on
  the unlogged intensity scale, so 4 is very small).

  A script for doing this comparison may be found in the \code{scripts}
  directory (although it is of limited use without the data files).
}
\value{
  \code{preprocessIllumina} returns a \code{MethylSet}, while
  \code{bgcorrect.illumina} and \code{normalize.illumina.control} both
  return a \code{RGChannelSet} with corrected color channels.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\linkS4class{RGChannelSet}} and \code{\linkS4class{MethylSet}}
  as well as \code{\linkS4class{IlluminaMethylationManifest}} for the
  basic classes involved in these functions.
  \code{\link{preprocessRaw}} is another basic preprocessing function.
}
\examples{
if (require(minfiData)) {

dat <- preprocessIllumina(RGsetEx, bg.correct=FALSE, normalize="controls")
slot(name="preprocessMethod", dat)[1]

}

}
