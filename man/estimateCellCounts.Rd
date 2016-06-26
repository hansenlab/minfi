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
estimateCellCounts(rgSet, compositeCellType = "Blood",
                   processMethod = "auto", probeSelect = "auto",
                   cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"),
                   returnAll = FALSE, meanPlot = FALSE, verbose = TRUE, \dots)
}
\arguments{
  \item{rgSet}{The input \code{RGChannelSet} for the procedure.}
  \item{compositeCellType}{Which composite cell type is being deconvoluted. Should
	be one of "Blood", "CordBlood", or "DLPFC". See details. }
  \item{processMethod}{How should the user and reference data be processed together?
	Default input "auto" will use \code{preprocessQuantile} for
	Blood and DLPFC and \code{preprocessNoob} otherwise, in line
	with the existing literature.  Set it to the name of a
	preprocessing function as a character if you want to override
	it, like \code{"preprocessFunnorm"}.}
  \item{probeSelect}{How should probes be selected to distinguish cell types? Options include 
	"both", which selects an equal number (50) of probes (with F-stat p-value < 1E-8) with the 
	greatest magnitude of effect from the hyper- and hypo-methylated sides, and "any", which 
	selects the 100 probes (with F-stat p-value < 1E-8) with the greatest magnitude of difference 
	regardless of direction of effect. Default input "auto" will use "any" for cord blood and 
	"both" otherwise, in line with previous versions of this function and/or our recommendations. 
	Please see the references for more details.}
  \item{cellTypes}{Which cell types, from the reference object, should be
    we use for the deconvolution? See details. }
  \item{returnAll}{Should the composition table and the normalized user
    supplied data be return?}
  \item{verbose}{Should the function be verbose?}
  \item{meanPlot}{
	Whether to plots the average DNA methylation across the cell-type discrimating probes within
		the mixed and sorted samples. 
  }
  \item{\dots}{Passed to \code{preprocessQuantile}.}
}
\details{

This is an implementaion of the Houseman et al (2012) regression calibration approachalgorithm
to the Illumina 450k microarray for deconvoluting heterogeneous tissue sources like blood.
For example, this function will take an \code{RGChannelSet} from a DNA methylation (DNAm)
study of blood, and return the relative proportions of CD4+ and CD8+ T-cells, natural
killer cells, monocytes,  granulocytes, and b-cells in each sample.

The function currently supports cell composition estimation for blood, cord blood, and 
the frontal cortex, through \code{compositeCellType} values of "Blood", "CordBlood", and 
"DLPFC", respectively. Packages containing the appropriate reference data should be installed 
before running the function for the first time ("FlowSorted.Blood.450k", "FlowSorted.DLPFC.450k",
"FlowSorted.CordBlood.450k"). Each tissue supports the estimation of different cell types, delimited 
via the \code{cellTypes} argument. For blood, these are "Bcell", "CD4T", "CD8T", "Eos", "Gran",
"Mono", "Neu", and "NK" (though the default value for \code{cellTypes} is often sufficient).
For cord blood, these are "Bcell", "CD4T", "CD8T", "Gran", "Mono", "Neu", and "nRBC". For frontal 
cortex, these are "NeuN_neg" and "NeuN_pos". See documentation of individual reference packages for
more details. 

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
  Matrix of composition estimates across all samples and cell types.

  If \code{returnAll=TRUE} a list of a count matrix (see previous
  paragraph), a composition table and the normalized user data in form
  of a GenomicMethylSet.
}

\references{
  EA Houseman, WP Accomando, DC Koestler, BC Christensen, CJ Marsit, HH
  Nelson, JK Wiencke and KT Kelsey.
  \emph{DNA methylation arrays as surrogate measures of cell mixture distribution.}
  BMC bioinformatics (2012) 13:86.
  doi:\href{http://www.dx.doi.org10.1186/1471-2105-13-86}{10.1186/1471-2105-13-86}.
  
  AE Jaffe and RA Irizarry. \emph{Accounting for cellular
    heterogeneity is critical in epigenome-wide association studies.}
  Genome Biology (2014) 15:R31.
  doi:\href{http://www.dx.doi.org/10.1186/gb-2014-15-2-r31}{10.1186/gb-2014-15-2-r31}.
  
  KM Bakulski, JI Feinberg, SV Andrews, J Yang, S Brown, S McKenney, 
  F Witter, J Walston, AP Feinberg, and MD Fallin. 
  \emph{DNA methylation of cord blood cell types: Applications for mixed cell birth studies.}
  Manuscript in review.
}

\author{
Andrew E. Jaffe, Shan V. Andrews, E. Andres Houseman
}

\seealso{\code{\link{preprocessQuantile}}}

\examples{
\dontrun{
if(require(FlowSorted.Blood.450k)) {
  wh.WBC <- which(FlowSorted.Blood.450k$CellType == "WBC")
  wh.PBMC <- which(FlowSorted.Blood.450k$CellType == "PBMC")
  RGset <- FlowSorted.Blood.450k[, c(wh.WBC, wh.PBMC)]
  ## The following line is purely to work around an issue with repeated
  ## sampleNames and Biobase::combine()
  sampleNames(RGset) <- paste(RGset$CellType,
    c(seq(along = wh.WBC), seq(along = wh.PBMC)), sep = "_")
  counts <- estimateCellCounts(RGset, meanPlot = FALSE)
  round(counts, 2)
}
}
}