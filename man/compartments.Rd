\name{compartments}
\alias{compartments}
\alias{createCorMatrix}
\alias{extractAB}
\title{
  Estimates A/B compartments from Illumina methylation arrays
}
\description{
  Estimates A/B compartments as revealed by Hi-C by computing the first
  eigenvector on a binned probe correlation matrix.
}
\usage{
compartments(object, resolution=100*1000, what = "OpenSea", chr="chr22",
                  method = c("pearson", "spearman"), keep=TRUE)
}
\arguments{
  \item{object}{An object of class (Genomic)MethylSet or (Genomic)RatioSet}
  \item{resolution}{An integer specifying the binning resolution}
  \item{what}{Which subset of probes should be used?}
  \item{chr}{The chromosome to be analyzed.}
  \item{method}{Method of correlation.}
  \item{keep}{Should the correlation matrix be stored or not?}
}
\details{
  
This function extracts A/B compartments from Illumina methylation microarrays. Analysis of Hi-C data has shown that the genome can be divided into two compartments (A/B compartments) that are cell-type specific and are associated with open and closed chromatin respectively. The approximately 170,000 open sea probes on the 450k array can be used to estimate these compartments by computing the first eigenvector on a binned correlation matrix. The binning resolution can be specified by \code{resolution}, and by default is set to a 100 kb. We do not recommend higher resolutions because of the low-resolution probe design of the 450k array.  
}

\value{an object of class \code{GRanges} containing the correlation matrix, the compartment eigenvector and the compartment labels (A or B) as metadata.}

\references{
  JP Fortin and KD Hansen.  
  \emph{Reconstructing A/B compartments as revealed by Hi-C using long-range correlations in epigenetic data}.
  bioRxiv (2015).
  doi:\href{http://dx.doi.org/10.1101/019000}{10.1101/019000}.
}

\author{
  Jean-Philippe Fortin \email{jfortin@jhsph.edu},
  Kasper D. Hansen \email{kasperdanielhansen@gmail.com}
}
\examples{
if (require(minfiData)) {
  GMset <- mapToGenome(MsetEx)
  comps <- compartments(GMset)
}
}
