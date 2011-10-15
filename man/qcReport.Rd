\name{qcReport}
\alias{qcReport}
\title{
  QC report for Illumina Infinium Human Methylation 450k arrays 
}
\description{
  Produces a PDF QC report for Illumina Infinium Human Methylation 450k
  arrays, useful for identifying failed samples. 
}
\usage{
qcReport(rgSet, sampNames = NULL, sampGroups = NULL, pdf = "qcReport.pdf",
    maxSamplesPerPage = 24, controls = c("BISULFITE CONVERSION I",
    "BISULFITE CONVERSION II", "EXTENSION", "HYBRIDIZATION",
    "NON-POLYMORPHIC", "SPECIFICITY I", "SPECIFICITY II", "TARGET REMOVAL"))
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet}.}
  \item{sampNames}{Sample names to be used for labels.}
  \item{sampGroups}{Sample groups to be used for labels.}
  \item{pdf}{Path and name of the PDF output file.}
  \item{maxSamplesPerPage}{Maximum number of samples to plot per page in
    those sections that plot each sample separately.} 
  \item{controls}{The control probe types to include in the report.}
}
\details{

  This function produces a QC report as a PDF file. It is a useful first
  step after reading in a new dataset to get an overview of quality and
  to flag potentially problematic samples.

}
\value{
No return value. A PDF is produced as a side-effect.
}
\author{
Martin Aryee \email{aryee@jhu.edu}.
}

\seealso{
\code{\link{mdsPlot}}, \code{\link{controlStripPlot}}, \code{\link{densityPlot}}, \code{\link{densityBeanPlot}}
}
\examples{
if (require(minfiData)) {

names <- pData(RGsetEx)$Sample_Name
groups <- pData(RGsetEx)$Sample_Group

\dontrun{
qcReport(RGsetEx, sampNames=names, sampGroups=groups, pdf="qcReport.pdf")
}

}
}