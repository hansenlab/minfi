\name{noobPipeline}
\alias{noobPipeline}
\alias{noobPipeline}
\title{
  A TCGA-style preprocessing pipeline for Infinium methylation microarrays.
}
\description{
  The standard processing pipeline for the Cancer Genome Atlas (TCGA) was 
  the methylumi implementation of preprocessNoob (which returns identical beta/M
  values to the ssNoob implementation now standard in minfi).  Probes with 
  detection p-values > 0.05 were masked as NA, SNP probes for label swap QC
  tabulated, and the results uploaded as Level 2 data after certain masking 
  steps (related to common SNPs and repetitive regions) had been applied.
  This function wraps the individual steps to turn an \code{RGChannelSet} into 
  a \code{GenomicRatioSet} directly comparable to canonical TCGA datasets. 
}
\usage{
noobPipeline(rgSet, pCutoff=0.05)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet}.}
  \item{pCutoff}{Above what detection pvalue should probes be marked NA? (0.05)}
}
\value{
  An object of class \code{GenomicRatioSet} with additional SNP metadata.
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
}
\examples{
if (require(minfiData)) {
  grSetEx <- noobPipeline(RGsetEx)
  identical(getSNPs(grSetEx), getSNPs(RGsetEx))
}
}
