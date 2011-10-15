\name{read.450k.sheet}
\alias{read.450k.sheet}
\title{
  Reading an Illumina methylation sample sheet
}
\description{
  Reading an Illumina methylation sample sheet, containing pheno-data
  information for the samples in an experiment.
}
\usage{
read.450k.sheet(base, pattern = "csv$", ignore.case = TRUE,
    recursive = TRUE, verbose = TRUE)
}
\arguments{
  \item{base}{The base directory from which the search is started.}
  \item{pattern}{What pattern is used to identify a sample sheet file,
    see \code{list.files}}
  \item{ignore.case}{Should the file search be case sensitive?}
  \item{recursive}{Should the file search be recursive, see \code{list.files}?}
  \item{verbose}{Should the function be verbose?}
}
\details{
    This function search the directory \code{base} (possibly including
  subdirectories depending on the argument \code{recursive} for
  \dQuote{sample sheet} files (see below).  These files are identified
  solely on the base of their filename given by the arguments
  \code{pattern} and \code{ignore.case} (note the use of a dollarsign to
  mean end of file name).

  In case multiple sheet files are found, they are all read and the
  return object will contain the concatenation of the files.
  
  A sample sheet file is essentially a CSV (comma-separated) file
  containing one line per sample, with a number of columns describing
  pheno-data or other important information about the sample.  The file
  may contain a header, in which case it is assumed that all lines up to
  and including a line starting with \code{\[Data\]} should be dropped.
  This is modelled after a sample sheet file Illumina provides.  It is
  also very similar to the \code{targets} file made used by the popular
  \code{limma} package (see the extensive package vignette).

  An attempt at guessing the file path to the IDAT files represented in
  the sheet is made.  This should be doublechecked and might need to
  manually changed.
}
\value{
  A \code{data.frame} containing the columns of all the sample sheets.
  As described in details, a column named \code{Sentrix_Position} is
  renamed to \code{Array} and \code{Sentrix_ID} is renamed to
  \code{Slide}.  In addition the \code{data.frame} will contain a column
  named \code{Basename}.
}
\author{
  Kasper Daniel Hansen\email{khansen@jhsph.edu}.
}
\seealso{
  \code{\link{read.450k.exp}} and \code{\link{read.450k}} for functions
  reading IDAT files.  \code{\link[base]{list.files}} for help on the
  arguments \code{recursive} and \code{ignore.case}.
}
\examples{
if(require(minfiData)) {

baseDir <- system.file("extdata", package = "minfiData")
sheet <- read.450k.sheet(baseDir)

}
}
