\name{densityPlot}
\alias{densityPlot}
\title{
Density plots of methylation Beta values.
}
\description{
Density plots of methylation Beta values, primarily for QC.
}
\usage{
densityPlot(dat, sampGroups = NULL, main = "", xlab = "Beta",
    pal = brewer.pal(8, "Dark2"), xlim, ylim, add = TRUE, legend = TRUE,
    \dots)
}
\arguments{
  \item{dat}{An \code{RGChannelSet}, a \code{MethylSet} or a
    \code{matrix}.  We either use the \code{getBeta} function to get
    Beta values (for the first two) or we assume the matrix contains
    Beta values.}
  \item{sampGroups}{Optional sample group labels. See details.}
  \item{main}{Plot title.}
  \item{xlab}{x-axis label.}
  \item{pal}{Color palette.}
  \item{xlim}{x-axis limits.}
  \item{ylim}{y-axis limits.}
  \item{add}{Start a new plot?}
  \item{legend}{Plot legend.}
  \item{\dots}{Additional options to be passed to the plot command.}
}
\details{
This function produces the density plot component of the QC report. If \code{sampGroups} is specified, group-specific colors will be used.
}
\value{
No return value. Plots are produced as a side-effect.
}
\author{
Martin Aryee \email{aryee@jhu.edu}.
}
\seealso{
\code{\link{qcReport}}, \code{\link{mdsPlot}}, \code{\link{controlStripPlot}}, \code{\link{densityBeanPlot}}
}
\examples{
if (require(minfiData)) {

groups <- pData(RGsetEx)$Sample_Group
densityPlot(RGsetEx, sampGroups=groups)

}
}
