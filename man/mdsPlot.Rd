\name{mdsPlot}
\alias{mdsPlot}
\title{
Multi-dimensional scaling plots giving an overview of similarities and differences between samples.
}
\description{
Multi-dimensional scaling (MDS) plots showing a 2-d projection of distances between samples.
}
\usage{
mdsPlot(dat, numPositions = 1000, sampNames = NULL, sampGroups = NULL, xlim, ylim,
    pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
    legendNCol, main = NULL)
}
\arguments{
  \item{dat}{An \code{RGChannelSet}, a \code{MethylSet} or a
    \code{matrix}.  We either use the \code{getBeta} function to get
    Beta values (for the first two) or we assume the matrix contains
    Beta values.}
  \item{numPositions}{Use the \code{numPositions} genomic positions with the most methylation variability when calculating distance between samples.}
  \item{sampNames}{Optional sample names. See details.}
  \item{sampGroups}{Optional sample group labels. See details.}
  \item{xlim}{x-axis limits.}
  \item{ylim}{y-axis limits.}
  \item{pch}{Point type. See \code{\link[graphics]{par}} for details.}
  \item{pal}{Color palette.}
  \item{legendPos}{The legend position. See
    \code{\link[graphics]{legend}} for details.} 
  \item{legendNCol}{The number of columns in the legend. See
    \code{\link[graphics]{legend}} for details.} 
  \item{main}{Plot title.}
}
\details{

  Euclidean distance is calculated between samples using the
  \code{numPositions} most variable CpG positions. These distances are then
  projected into a 2-d plane using classical multidimensional scaling
  transformation.

}
\value{
No return value. Plots are produced as a side-effect.
}
\references{

  Borg, I., Groenen, P. (2005). Modern Multidimensional Scaling: theory
and applications (2nd ed.). New York: Springer-Verlag. pp. 207-212. ISBN
0387948457.

  \url{http://en.wikipedia.org/wiki/Multidimensional_scaling}
}
\author{
Martin Aryee \email{aryee@jhu.edu}.
}

\seealso{
\code{\link{qcReport}}, \code{\link{controlStripPlot}},
\code{\link{densityPlot}}, \code{\link{densityBeanPlot}},
\code{\link[graphics]{par}}, \code{\link[graphics]{legend}}
}
\examples{
if (require(minfiData)) {

names <- pData(MsetEx)$Sample_Name
groups <- pData(MsetEx)$Sample_Group
mdsPlot(MsetEx, sampNames=names, sampGroups=groups)

}
}