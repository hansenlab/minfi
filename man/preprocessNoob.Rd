\name{preprocessNoob}
\alias{preprocessNoob}
\title{
  The Noob/ssNoob preprocessing method for Infinium methylation microarrays.
}
\description{
  Noob (normal-exponential out-of-band) is a background correction method
  with dye-bias normalization for Illumina Infinium methylation arrays.
}
\usage{
preprocessNoob(rgSet, offset = 15, dyeCorr = TRUE, verbose = FALSE,
               dyeMethod=c("single", "reference"))
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet}.}
  \item{offset}{An offset for the normexp background correction.}
  \item{dyeCorr}{Should dye correction be done?}
  \item{verbose}{Should the function be verbose?}
  \item{dyeMethod}{How should dye bias correction be done: use a single sample approach (ssNoob), or a reference array?}
}
\value{
  An object of class \code{MethylSet}.
}
\references{
  TJ Triche, DJ Weisenberger, D Van Den Berg, PW Laird and KD Siegmund
  \emph{Low-level processing of Illumina Infinium DNA Methylation
    BeadArrays}.
  Nucleic Acids Res (2013) 41, e90.
  doi:\href{http://www.dx.doi.org/10.1093/nar/gkt090}{10.1093/nar/gkt090}.
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
  ## RGsetEx.sub is a small subset of RGsetEx;
  ## only used for computational speed.
  MsetEx.sub.noob <- preprocessNoob(RGsetEx.sub)
}
\dontrun{
if (require(minfiData)) {
  dyeMethods <- c(ssNoob="single", refNoob="reference")
  GRsets <- lapply(dyeMethods,
                   function(m) preprocessNoob(RGsetEx, dyeMethod=m))
  all.equal(getBeta(GRsets$refNoob), getBeta(GRsets$ssNoob)) # TRUE
}
}
}
