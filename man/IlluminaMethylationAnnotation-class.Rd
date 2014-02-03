\name{IlluminaMethylationAnnotation-class}
\Rdversion{1.1}
\docType{class}
\alias{IlluminaMethylationAnnotation-class}
\alias{IlluminaMethylationAnnotation}
\alias{show,IlluminaMethylationAnnotation-method}
\alias{getManifest,IlluminaMethylationAnnotation-method}

\title{Class \code{IlluminaMethylationAnnotation}}

\description{
  This is a class for representing annotation associated with an
  Illumina methylation microarray.  Annotation is transient in the sense
  that it may change over time, wheres the information stored in the
  \code{IlluminaMethylationManifest} class only depends on the array design.
}

\usage{
## Constructor
IlluminaMethylationAnnotation(listOfObjects, annotation = "",
                              defaults = "")

## Data extraction
\S4method{getManifest}{IlluminaMethylationAnnotation}(object)
}

\arguments{
  \item{object}{An object of class
    \code{IlluminaMethylationAnnotation}.}
  \item{annotation}{An annotation \code{character}.}
  \item{defaults}{A vector of default choices for
    \code{getAnnotation(what = "everything")}.}
  \item{listOfObjects}{A list of objects to be put into the data slot of
    the annotation object.}
}

\section{Details}{
  The \code{data} slot contains various objects that are
  of type \code{data.frame}.  Details are still subject to change. 
}

\section{Utilitues}{
  In the following code, \code{object} is a \code{IlluminaMethylationAnnotation}.

  \describe{
    \item{\code{getManifest(object)}}{Get the manifest object associated
      with the array.}
  }
}

\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\linkS4class{IlluminaMethylationManifest}}
}
\keyword{classes}
