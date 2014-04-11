\name{dmpFinder}
\alias{dmpFinder}
\title{
Find differentially methylated positions
}
\description{
Identify CpGs where methylation is associated with a continuous or categorical phenotype.
}
\usage{
dmpFinder(dat, pheno, type = c("categorical", "continuous"),
    qCutoff = 1, shrinkVar = FALSE)
}
\arguments{
  \item{dat}{ A \code{MethylSet} or a \code{matrix}.}
  \item{pheno}{ The phenotype to be tested for association with methylation. }
  \item{type}{ Is the phenotype '\sQuote{continuous} or \sQuote{categorical}?  }
  \item{qCutoff}{ DMPs with an FDR q-value greater than this will not be returned.  }
  \item{shrinkVar}{ Should variance shrinkage be used? See details. }
}
\details{

  This function tests each genomic position for association between methylation and
  a phenotype. Continuous phenotypes are tested with linear regression,
  while an F-test is used for categorical phenotypes. 

  Variance shrinkage (\code{shrinkVar=TRUE}) is recommended when sample
  sizes are small (<10). The sample variances are squeezed
  by computing empirical Bayes posterior means using the \pkg{limma}
  package.

}
\value{
  A table with one row per CpG.
}
\author{
  Martin Aryee \email{aryee@jhu.edu}.
}
\seealso{
  \code{\link[limma]{squeezeVar}} and the \pkg{limma} package in general.
}
\examples{
if (require(minfiData)) {

grp <- pData(MsetEx)$Sample_Group
MsetExSmall <- MsetEx[1:1e4,] # To speed up the example
M <- getM(MsetExSmall, type = "beta", betaThreshold = 0.001)
dmp <- dmpFinder(M, pheno=grp, type="categorical")
sum(dmp$qval < 0.05, na.rm=TRUE)
head(dmp)

}
}