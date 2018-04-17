\name{getSex}
\alias{getSex}
\alias{addSex}
\alias{plotSex}
\title{
  Estimating sample sex based on methylation data
}
\description{
  Estimates samples sex based on methylation data.
}
\usage{
getSex(object = NULL, cutoff = -2)
addSex(object, sex = NULL)
plotSex(object, id = NULL)
}
\arguments{
  \item{object}{An object of class \code{[Genomic]MethylSet}.}
  \item{cutoff}{What should the difference in log2 copynumber be between
  males and females.}
  \item{sex}{An optional character vector of sex (with values \code{M}
    and \code{F}).}
  \item{id}{Text used as plotting symbols in the \code{plotSex}
    function. Used for sample identification on the plot.}
}
\details{
  Estimation of sex is based on the median values of measurements on the
  X and Y chromosomes respectively. If \code{yMed} - \code{xMed} is
  less than \code{cutoff} we predict a female, otherwise male.
}
\value{
  For \code{getSex}, a \code{DataFrame} with columns \code{predictedSex}
  (a character with values \code{M} and \code{F}), \code{xMed} and
  \code{yMed}, which are the chip-wide medians of measurements on the two
  sex chromosomes.

  For \code{addSex}, an object of the same type as \code{object} but
  with the output of \code{getSex(object)} added to the pheno data.

  For \code{plotSex}, a plot of \code{xMed} vs. \code{yMed}, which are the
  chip-wide medians of measurements on the two sex chromosomes, coloured by
  \code{predictedSex}.

}
\author{
  Rafael A. Irizarry, Kasper D. Hansen, Peter F. Hickey
}
\examples{
if(require(minfiData)) {
  GMsetEx <- mapToGenome(MsetEx)
  estSex <- getSex(GMsetEx)
  GMsetEx <- addSex(GMsetEx, sex = estSex)
}
}
