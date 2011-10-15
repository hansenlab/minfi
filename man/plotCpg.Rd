\name{plotCpg}
\alias{plotCpg}
\title{
Plot methylation values at an single genomic position
}
\description{
Plot single-position (single CpG) methylation values as a function of a categorical or continuous phenotype
}
\usage{
plotCpg(dat, cpg, pheno, type = c("categorical", "continuous"),
    measure = c("beta", "M"), ylim = NULL, ylab = NULL, xlab = "",
    fitLine = TRUE, mainPrefix = NULL, mainSuffix = NULL) 
}
\arguments{
    \item{dat}{An \code{RGChannelSet}, a \code{MethylSet} or a
    \code{matrix}.  We either use the \code{getBeta} (or \code{getM} for
    \code{measure="M"}) function to get Beta values (or M-values) (for
    the first two) or we assume the matrix contains Beta values (or
    M-values).}
  \item{cpg}{A character vector of the CpG position identifiers to be plotted.}
  \item{pheno}{A vector of phenotype values.}
  \item{type}{Is the phenotype categorical or continuous?}
  \item{measure}{Should Beta values or log-ratios (M) be  plotted?}
  \item{ylim}{y-axis limits.}
  \item{ylab}{y-axis label.}
  \item{xlab}{x-axis label.}
  \item{fitLine}{Fit a least-squares best fit line when using a continuous phenotype.}
  \item{mainPrefix}{Text to prepend to the CpG name in the plot main title.}
  \item{mainSuffix}{Text to append to the CpG name in the plot main title.}
}
\details{
  This function plots methylation values (Betas or log-ratios) at
  individual CpG loci as a function of a phenotype.

}
\value{
No return value. Plots are produced as a side-effect.
}

\author{
Martin Aryee \email{aryee@jhu.edu}.
}

\examples{
if (require(minfiData)) {

grp <- pData(MsetEx)$Sample_Group
cpgs <- c("cg00050873", "cg00212031", "cg26684946", "cg00128718")
par(mfrow=c(2,2))
plotCpg(MsetEx, cpg=cpgs, pheno=grp, type="categorical")

}
}

