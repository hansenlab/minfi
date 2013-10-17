\name{estimateCellCounts}
\alias{estimateCellCounts}
\title{
Cell Proportion Estimation
}
\description{
Estimates the relative proportion of pure cell types within a sample. For example, 
	given peripheral blood samples, this function will return the relative proportions
	of lymphocytes, monocytes, B-cells, and neutrophils. 
}
\usage{
estimateCellCounts(RGset, meanPlot = TRUE, cellType = "Blood", verbose = TRUE, \dots)
}
\arguments{
  \item{RGset}{The input \code{RGChannelSet} for the procedure.}
  \item{meanPlot}{
	Whether to plots the average DNA methylation across the cell-type discrimating probes within
		the mixed and sorted samples. 
  }
  \item{cellType}{Which cell type is being deconvoluted, see details.}
  \item{verbose}{Should the function be verbose?}
  \item{\dots}{Passed to \code{preprocessQuantile}.}
}
\details{

This is an implementaion of the Houseman et al (2012) regression calibration approachalgorithm
to the Illumina 450k microarray for deconvoluting heterogeneous tissue sources like blood.
For example, this function will take an \code{RGChannelSet} from a DNA methylation (DNAm)
study of blood, and return the relative proportions of CD4+ and CD8+ T-cells, natural
killer cells, monocytes,  granulocytes, and b-cells in each sample.

The \code{meanPlot} should be used to check for large batch effects in the data,
reducing the confidence placed in the composition estimates. This plot
depicts the average DNA methylation across the cell-type discrimating probes
in both the provided and sorted data. The means from the provided 
heterogeneous samples should be within the range of the sorted samples.
If the sample means fall outside the range of the sorted means, 
the cell type estimates will inflated to the closest cell type. Note that we
quantile normalize the sorted data with the provided data to reduce these
batch effects. 

}
\value{
  Matrix of composition estimates across all samples and cell types
}

\references{
Houseman, Eugene A., et al. "DNA methylation arrays as surrogate
measures of cell mixture distribution." BMC bioinformatics 13.1 (2012):
86. 
	
Andrew E. Jaffe and Rafael A. Irizarry. Accounting for cellular
heterogeneity is critical in epigenome-wide association studies. Under
Review.
}

\author{
Andrew E. Jaffe, E. Andres Houseman
}

\note{
  At the time of writing, only data from deconvolution of whole blood
  cells are available.
}

\seealso{\code{\link{preprocessQuantile}}}
