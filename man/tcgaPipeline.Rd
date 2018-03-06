\name{tcgaPipeline}
\alias{tcgaPipeline}
\alias{tcgaPipeline}
\title{
  A TCGA-style preprocessing pipeline for Infinium methylation microarrays.
}
\description{
  The standard DNAme processing pipeline for the Cancer Genome Atlas (TCGA) was 
  the methylumi implementation of preprocessNoob (which yields identical beta/M
  values to the single-sample ssNoob implementation now standard in minfi, in 
  the common case where a balanced reference sample was identified in each TCGA
  batch). Probes with detection p-values > 0.05 were masked as NA, and SNP (rs)
  probes included on each array by Illumina for label swap detection tabulated.
  After QC and certain masking steps (related to common SNPs and repetitive 
  regions) had been applied, the results were uploaded as Level 2 data for TCGA.
  This function wraps the individual steps to turn an \code{RGChannelSet} into 
  a \code{GenomicRatioSet} directly comparable to canonical TCGA datasets, save 
  for the SNP/repeat masking steps (which, in hindsight, seem overzealous). The 
  resulting GenomicRatioSet is thus directly comparable to processed TCGA data,
  and as of this writing (May 2017), no single preprocessing method (let alone 
  any single-sample preprocessing method) has emerged as superior for the
  general case of mixed tumor and normal samples from varying tissues. Thus, 
  this approach may also serve as a useful benchmark for methods development.
}
\usage{
tcgaPipeline(rgSet, pCutoff=0.05)
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
  grSetEx <- tcgaPipeline(RGsetEx)
  stopifnot(identical(getSnpBeta(grSetEx), getSnpBeta(RGsetEx)))
}
}
