\name{IlluminaMethylationAnnotation-class}
\Rdversion{1.1}
\docType{class}
\alias{IlluminaMethylationAnnotation-class}
\alias{show,IlluminaMethylationAnnotation-method}
\title{Class \code{"IlluminaMethylationAnnotation"}}
\description{
  This is a class for representing annotation associated with an
  Illumina methylation microarray.  Annotation is transient in the sense
  that it may change over time, wheres the information stored in the
  \code{"IlluminaMethylationManifest"} class only depends on the array design.
}
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{new("IlluminaMethylationAnnotation", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{This environment holds the annotation objects.  See details.}
    \item{\code{annotation}:}{A character vector giving the annotation name.}
  }
}
\section{Details}{
  The \code{data} slot contains various objects that are
  of type \code{data.frame}.  Details are still subject to change. 
}
\section{Methods}{
  \describe{
    \item{\code{show}:}{The show method.}
  }
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\linkS4class{IlluminaMethylationManifest}}
}
\keyword{classes}
