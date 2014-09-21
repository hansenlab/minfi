\name{preprocessNoob}
\alias{preprocessNoob}
\title{
  The Noob background correction for 450k arrays.
}
\description{
  Noob is a background correction method with dye-bias normalization for
  the Illumina Infinium HumanMethylation450 platform.
}
\usage{
preprocessNoob(rgSet, offset = 15, dyeCorr = TRUE, verbose = TRUE)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet}.}
  \item{offset}{An offset for the normexp background correction.}
  \item{dyeCorr}{Should dye correction be done?}
  \item{verbose}{Should the function be verbose?}
}
\value{
  An object of class \code{MethylSet}.
}
\references{
  TJ Triche, DJ Weisenberger, D Van Den Berg, PW Laird and KD Siegmund
  \emph{Low-level processing of Illumina Infinium DNA Methylation
    BeadArrays}.
  (2013) Nucleic Acids Res 41, e90.
}
\author{
  Tim Triche, Jr.
}
\seealso{
  \code{\linkS4class{RGChannelSet}}
  as well as \code{\linkS4class{IlluminaMethylationManifest}} for the
  basic classes involved in these functions.
  \code{\link{preprocessRaw}} and \code{\link{preprocessQuantile}} are other preprocessing functions.
}
\examples{
if (require(minfiData)) {
  GRset.noob <- preprocessNoob(RGsetEx)
}
}
