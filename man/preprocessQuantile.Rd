\name{preprocessQuantile}
\alias{preprocessQuantile}
\title{
  Stratified quantile normalization for an Illumina methylation array.
}

\description{
  Stratified quantile normalization for Illumina amethylation arrays.
  
  This function implements stratified quantile normalization preprocessing for 
  Illumina methylation microarrays. Probes are stratified by region (CpG island,
  shore, etc.)

}
\usage{
preprocessQuantile(object, fixOutliers = TRUE, removeBadSamples = FALSE,
                   badSampleCutoff = 10.5, quantileNormalize = TRUE,
                   stratified = TRUE, mergeManifest = FALSE, sex = NULL,
                   verbose = TRUE)
}
\arguments{
  \item{object}{An object of class \code{RGChannelSet} or
    \code{[Genomic]MethylSet}.}
  \item{fixOutliers}{Should low outlier Meth and Unmeth signals be
    fixed?} 
  \item{removeBadSamples}{Should bad samples be removed?}
  \item{badSampleCutoff}{Samples with median Meth and Umneth signals
    below this cutoff will be labelled \sQuote{bad}.} 
  \item{quantileNormalize}{Should quantile normalization be performed?}
  \item{stratified}{Should quantile normalization be performed within
    genomic region strata (e.g. CpG island, shore, etc.)?}
  \item{mergeManifest}{Should the information in the associated manifest
    package be merged into the output object?}
  \item{sex}{Gender}
  \item{verbose}{Should the function be verbose?}
}
\details{

This function implements stratified quantile normalization preprocessing
for Illumina methylation microarrays. If \code{removeBadSamples} is
\code{TRUE} we calculate the median Meth and median Unmeth signal for
each sample, and remove those samples where their average falls below
\code{badSampleCutoff}. The normalization procedure is applied to the
Meth and Unmeth intensities separately. The distribution of type I and
type II signals is forced to be the same by first quantile normalizing
the type II probes across samples and then interpolating a reference
distribution to which we normalize the type I probes. Since probe types
and probe regions are confounded and we know that DNAm distributions
vary across regions we stratify the probes by region before applying
this interpolation. For the probes on the X and Y chromosomes we
normalize males and females separately using the gender information
provided in the \code{sex} argument. If gender is unspecified
(\code{NULL}), a guess is made using by the \code{getSex} function using
copy number information. Background correction is not used, but very
small intensities close to zero are thresholded using the
\code{fixMethOutlier}. Note that this algorithm relies on the
assumptions necessary for quantile normalization to be applicable and
thus is not recommended for cases where global changes are expected such
as in cancer-normal comparisons.

Note that this normalization procedure is essentially similar to one
previously presented (Touleimat and Tost, 2012), but has been
independently re-implemented due to the present lack of a released,
supported version.

}
\value{
  a \code{GenomicRatioSet}
}

\note{
  A bug in the function was found to affect the Beta values of type I
  probes, when \code{stratified=TRUE} (default). This is fixed in minfi
  version 1.19.7 and 1.18.4 and greater.
}

\references{
  N Touleimat and J Tost.
  \emph{Complete pipeline for Infinium Human
    Methylation 450K BeadChip data processing using subset quantile
    normalization for accurate DNA methylation estimation.}
  Epigenomics (2012) 4:325-341.
}

\author{
  Rafael A. Irizarry
}

\seealso{
  \code{\link{getSex}}, \code{\link{minfiQC}},
  \code{\link{fixMethOutliers}} for functions used as part of
  \code{preprocessQuantile}.
  }

\examples{
\dontrun{
if(require(minfiData)) {
  GMset <- preprocessQuantile(RGsetEx)
}
}
}
