\name{densityBeanPlot}
\alias{densityBeanPlot}
\title{
  Density bean plots of methylation Beta values.
}
\description{
  Density \sQuote{bean} plots of methylation Beta values, primarily for QC.
}
\usage{
densityBeanPlot(dat, sampGroups = NULL, sampNames = NULL, main = NULL,
    pal = brewer.pal(8, "Dark2"), numPositions = 10000)
}
\arguments{
  \item{dat}{An \code{RGChannelSet}, a \code{MethylSet} or a
    \code{matrix}.  We either use the \code{getBeta} function to get
    Beta values (for the first two) or we assume the matrix contains
    Beta values.}
  \item{sampGroups}{Optional sample group labels. See details.}
  \item{sampNames}{Optional sample names. See details.}
  \item{main}{Plot title.}
  \item{pal}{Color palette.}
  \item{numPositions}{The density calculation uses \code{numPositions}
    randomly selected CpG positions.  If \code{NULL} use all positions.}
}
\details{
  This function produces the density bean plot component of the QC
  report. If \code{sampGroups} is specified, group-specific colors will
  be used. For speed reasons the plots are produced using a random
  subset of CpG positions. The number of positions used is specified by the
  \code{numPositions} option.
}
\value{
No return value. Plots are produced as a side-effect.

}
\references{
Kampstra, P. Beanplot: A boxplot alternative for visual comparison of distributions. Journal of Statistical Software 28, (2008). http://www.jstatsoft.org/v28/c01
}
\author{
Martin Aryee \email{aryee@jhu.edu}.
}

\seealso{
\code{\link{qcReport}}, \code{\link{mdsPlot}}, \code{\link{controlStripPlot}}, \code{\link{densityPlot}}
}
\examples{
if (require(minfiData)) {

names <- pData(RGsetEx)$Sample_Name
groups <- pData(RGsetEx)$Sample_Group
par(mar=c(5,6,4,2))
densityBeanPlot(RGsetEx, sampNames=names, sampGroups=groups)

}
}
