\name{bumphunter-methods}
\docType{methods}
\alias{bumphunter,GenomicRatioSet-method}

\title{Methods for function \code{bumphunter} in Package \pkg{minfi}}
\description{
 Methods for function \code{bumphunter} in package \pkg{minfi}.
}
\section{Methods}{
\describe{
\item{\code{signature(object = "GenomicRatioSet")}}{
  An object of class \code{GenomicRatioSet}
}
}}
\examples{
if(require(minfiData)) {
  gmSet <- preprocessQuantile(MsetEx)
  design <- model.matrix(~ gmSet$status)
  bumps <- bumphunter(gmSet, design = design, B = 0,
                      type = "Beta", cutoff = 0.25)
}
}

