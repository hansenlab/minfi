\name{getGenomicRatioSetFromGEO}
\alias{getGenomicRatioSetFromGEO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reading Illumina methylation array data from GEO.
}
\description{
Reading Illumina methylation array data from GEO.
}
\usage{
getGenomicRatioSetFromGEO(GSE = NULL, path = NULL, array = "IlluminaHumanMethylation450k", annotation = .default.450k.annotation, what = c("Beta", "M"), mergeManifest = FALSE, i = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GSE}{
The GSE ID of the dataset to be downloaded from GEO.
}
  \item{path}{
If data already downloaded, the path with soft files. Either \code{GSE} or \code{path}
are required.
}
  \item{array}{
Array name.
}
  \item{annotation}{
The feature annotation to be used. This includes the location of
features thus depends on genome build.
}
  \item{what}{
Are \code{Beta} or \code{M} values being downloaded.  
}
\item{mergeManifest}{
Should the Manifest be merged to the final object.
}
  \item{i}{
If the GEO download results in more than one dataset, it pickes entry \code{i}.
}
}
\details{
This function downloads data from GEO using
\code{\link[GEOquery]{getGEO}} from the \pkg{GEOquery} package. It then
returns a \linkS4class{GenomicRatioSet} object. Note that the rs probes
(used for genotyping) are dropped.
}
\value{
A \linkS4class{GenomicRatioSet} object.
}
\author{
  Tim Triche Jr. and 
  Rafael A. Irizarry\email{rafa@jimmy.harvard.edu}.
}
\seealso{
  If the data is already in memor you can use 
  \code{\link{makeGenomicRatioSetFromMatrix}}
}
\examples{
\dontrun{
mset=getGenomicRatioSetFromGEO("GSE42752")
}
}


