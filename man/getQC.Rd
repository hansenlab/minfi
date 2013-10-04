\name{getQC}
\alias{getQC}
\alias{addQC}
\alias{plotQC}
\title{
  Estimate sample-specific quality control (QC) for methylation data 
}
\description{
  Estimate sample-specific quality control (QC) for methylation data.
}
\usage{
getQC(object)
addQC(object, qc)
plotQC(qc, badSampleCutoff = 10.5)
}
\arguments{
  \item{object}{An object of class \code{[Genomic]MethylSet}.}
  \item{qc}{An object as produced by \code{getQC}.}
  \item{badSampleCutoff}{The cutoff for identifying a bad sample.}
}
\value{
  For \code{getQC}, a \code{DataFrame} with two columns: \code{mMed} and
  \code{uMed} which are the chipwide medians of the Meth and Unmeth
  channels.

  For \code{addQC}, essentially \code{object} supplied to the function,
  but with two new columns added to the pheno data slot: \code{uMed} and
  \code{mMed}.
}
\author{
  Rafael A. Irizarry and Kasper D. Hansen
}
\seealso{
  \code{\link{minfiQC}} for an all-in-one function.
}
\examples{
if(require(minfiData)){
  qc <- getQC(MsetEx)
  MsetEx <- addQC(MsetEx, qc = qc)
  ## plotQC(qc)
}
}
