\name{controlStripPlot}
\alias{controlStripPlot}
\title{
Plot control probe signals.
}
\description{
Strip plots are produced for each control probe type specified.
}
\usage{
controlStripPlot(rgSet, controls = c("BISULFITE CONVERSION I",
    "BISULFITE CONVERSION II"), sampNames = NULL, xlim = c(5, 17))
}
\arguments{
  \item{rgSet}{An \code{RGChannelSet}.}
  \item{controls}{A vector of control probe types to plot.}
  \item{sampNames}{Sample names to be used for labels.}
  \item{xlim}{x-axis limits.}
}
\details{
This function produces the control probe signal plot component of the QC report.
}
\value{
No return value. Plots are produced as a side-effect.
}
\author{
Martin Aryee \email{aryee@jhu.edu}.
}
\seealso{
\code{\link{qcReport}}, \code{\link{mdsPlot}}, \code{\link{densityPlot}}, \code{\link{densityBeanPlot}}
}
\examples{
if (require(minfiData)) {

names <- pData(RGsetEx)$Sample_Name
controlStripPlot(RGsetEx, controls=c("BISULFITE CONVERSION I"), sampNames=names)

}
}
