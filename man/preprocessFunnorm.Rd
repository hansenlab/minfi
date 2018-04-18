\name{preprocessFunnorm}
\alias{preprocessFunnorm}
\title{
  Functional normalization for Illumina 450k arrays
}
\description{
  Functional normalization (FunNorm) is a between-array normalization
  method for the Illumina Infinium HumanMethylation450 platform. It
  removes unwanted variation by regressing out variability explained by
  the control probes present on the array.
}
\usage{
preprocessFunnorm(rgSet, nPCs=2, sex = NULL, bgCorr = TRUE,
                  dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE,
                  verbose = TRUE)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet}.}
  \item{nPCs}{Number of principal components from the control probes PCA}
  \item{sex}{An optional numeric vector containing the sex of the
    samples.}
  \item{bgCorr}{Should the NOOB background correction be done, prior to
    functional normalization (see \code{preprocessNoob})}
  \item{dyeCorr}{Should dye normalization be done as part of the NOOB
    background correction (see \code{preprocessNoob})?}
  \item{keepCN}{Should copy number estimates be kept around? Setting to
    \code{FALSE} will decrease the size of the output object significantly.}
  \item{ratioConvert}{Should we run \code{ratioConvert}, ie. should the
    output be a \code{GenomicRatioSet} or should it be kept as a
    \code{GenomicMethylSet}; the latter is for experts.}
  \item{verbose}{Should the function be verbose?}
}
\details{

This function implements functional normalization preprocessing for
Illumina methylation microarrays. Functional normalization extends the
idea of quantile normalization by adjusting for known covariates
measuring unwanted variation. For the 450k array, the first k principal
components of the internal control probes matrix play the role of the
covariates adjusting for technical variation. The number k of principal
components can be set by the argument \code{nPCs}. By default
\code{nPCs} is set to 2, and have been shown to perform consistently
well across different datasets.  This parameter should only be modified
by expert users.  The normalization procedure is applied to the Meth and
Unmeth intensities separately, and to type I and type II signals
separately. For the probes on the X and Y chromosomes we normalize males
and females separately using the gender information provided in the
\code{sex} argument. For the Y chromosome, standard quantile
normalization is used due to the small number of probes, which results
in instability for functional normalization. If sex is unspecified
(\code{NULL}), a guess is made using by the \code{getSex} function using
copy number information. Note that this algorithm does not rely on any
assumption and therefore can be be applicable for cases where global
changes are expected such as in cancer-normal comparisons or tissue
differences.

}

\value{an object of class \code{GenomicRatioSet}, unless
  \code{ratioConvert=FALSE} in which case an object of class
  \code{GenomicMethylSet}.}

\references{
  JP Fortin, A Labbe, M Lemire, BW Zanke, TJ Hudson, EJ Fertig, CMT
  Greenwood and KD Hansen.
  \emph{Functional normalization of 450k methylation array data improves
    replication in large cancer studies}.
  (2014) Genome Biology (2014) 15:503.
  doi:\href{http://www.dx.doi.org/10.1186/s13059-014-0503-2}{10.1186/s13059-014-0503-2}.
}

\author{
  Jean-Philippe Fortin \email{jfortin@jhsph.edu},
  Kasper D. Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\linkS4class{RGChannelSet}}
  as well as \code{\linkS4class{IlluminaMethylationManifest}} for the
  basic classes involved in these functions.
  \code{\link{preprocessRaw}} and \code{\link{preprocessQuantile}} are
  other preprocessing functions. Background correction may be done using
  \code{\link{preprocessNoob}}.
}
\examples{
if (require(minfiData)) {
  ## RGsetEx.sub is a small subset of RGsetEx;
  ## only used for computational speed.
  Mset.sub.funnorm <- preprocessFunnorm(RGsetEx.sub)
}
}
