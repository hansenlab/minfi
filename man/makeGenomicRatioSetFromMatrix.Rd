\name{makeGenomicRatioSetFromMatrix}
\alias{makeGenomicRatioSetFromMatrix}
\title{
Make a GenomicRatioSet from a matrix
}
\description{
Make a GenomicRatioSet from a matrix.
}
\usage{
makeGenomicRatioSetFromMatrix(mat, rownames = NULL, pData = NULL, array = "IlluminaHumanMethylation450k", annotation = .default.450k.annotation, mergeManifest = FALSE, what = c("Beta", "M"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
The matrix that will be converted.
}
  \item{rownames}{
The feature IDs associated with the rows of \code{mat} that will be used
to match to the IlluminaHumanMethylation450k feature IDs.
}
  \item{pData}{
A \code{DataFrame} or \code{data.frame} describing the samples represented by the columns of
\code{mat}. If the rownames of the \code{pData} don't match the colnames
of \code{mat} these colnames will be changed. If \code{pData} is not
supplied, a minimal \code{DataFrame} is created.
}
  \item{array}{
Array name.
}
  \item{annotation}{
The feature annotation to be used. This includes the location of
features thus depends on genome build.
}
\item{mergeManifest}{
Should the Manifest be merged to the final object.
}
\item{what}{
Are \code{Beta} or \code{M} values being downloaded.  
}
}
\details{
Many 450K data is provided as csv files. This function permits you to
convert a matrix of values into the class that is used by functions such
as \code{bumphunter} and \code{blockFinder}. The rownmaes of \code{mat}
are used to match the 450K array features. Alternatively the rownames
can be supplied directly through \code{rownames}.
}
\value{
  A \linkS4class{GenomicRatioSet} object.
}
\author{
  Rafael A. Irizarry\email{rafa@jimmy.harvard.edu}.
}
\seealso{
  \code{\link{getGenomicRatioSetFromGEO}} is similar but reads data from GEO.
}
\examples{
mat <- matrix(10,5,2)
rownames(mat) <- c( "cg13869341", "cg14008030","cg12045430", "cg20826792","cg00381604")
grset <- makeGenomicRatioSetFromMatrix(mat)
}
