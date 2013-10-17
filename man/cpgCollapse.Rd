\name{cpgCollapse}
\alias{cpgCollapse}
\title{
  Collapse methylation values of adjacent CpGs into a summary value.
}

\description{

  This function groups adjacent loci into clusters with a specified
maximum gap between CpGs in the cluster, and a specified maximum cluster
width.  The loci within each cluster are summarized resulting in a
single methylation estimate per cluster.

}

\usage{
cpgCollapse(object, what = c("Beta", "M"), maxGap = 500,
            blockMaxGap = 2.5 * 10^5, maxClusterWidth = 1500,
            dataSummary = colMeans, na.rm = FALSE,
            returnBlockInfo = TRUE, verbose = TRUE, \dots)
}

\arguments{
  \item{object}{An object of class \code{[Genomic]MethylSet}.}
  \item{what}{Should operation be performed on the M-scale or
  Beta-scale?}
  \item{maxGap}{Maximum gap between CpGs in a cluster}
  \item{blockMaxGap}{Maximum block gap}
  \item{maxClusterWidth}{Maximum cluster width}
  \item{dataSummary}{Function used to summarize methylation across CpGs
    in the cluster.}
  \item{na.rm}{Should NAs be removed when summarizing? Passed on to the
    dataSummary function.}
  \item{returnBlockInfo}{Should the block annotation table be returned
    in addition to the block table?} 
  \item{verbose}{Should the function be verbose?}
  \item{\dots}{Passed on to getMethSignal and getCN. Can be used to
    specify } 
}
\details{

  This function is used as the first step of block-finding. It groups
  adjacent loci into clusters with a default maximum gap of 500bp and a
  maximum cluster width of 1,500bp.  The loci within each cluster are
  then summarized (using the mean by default) resulting in a single
  methylation estimate per cluster. Cluster estimates from open-sea
  probes are used in block-finding.

}

\value{

If returnBlockInfo is \code{FALSE}: a GenomicRatioSet of collapsed CpG clusters.

If returnBlockInfo is \code{TRUE}:
\item{object }{A GenomicRatioSet of collapsed CpG clusters}
\item{blockInfo }{A cluster annotation data frame}
}

\author{
Rafael Irizarry
}

\seealso{
  \code{\link{blockFinder}}
}

