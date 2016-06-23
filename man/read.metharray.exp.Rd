\name{read.metharray.exp}
\alias{read.metharray.exp}
\title{
  Reads an entire metharray experiment using a sample sheet
}
\description{
  Reads an entire methylation array experiment using a sample sheet or (optionally)
  a target like data.frame.
}
\usage{
read.metharray.exp(base = NULL, targets = NULL, extended = FALSE,
    recursive = FALSE, verbose = FALSE, force = FALSE)
}
\arguments{
  \item{base}{The base directory.}
  \item{targets}{A targets \code{data.frame}, see details}
  \item{extended}{Should the output of the function be a
    \code{"RGChannelSetExtended"} (default is \code{"RGChannelSet"}).}
  \item{recursive}{Should the search be recursive (see details)}
  \item{verbose}{Should the function be verbose?}
  \item{force}{Should reading different size IDAT files be forced? See
    the documentation for \code{\link{read.metharray}}}
}
\details{
  If the \code{targets} argument is \code{NULL}, the function finds all
  two-color IDAT files in the directory given by \code{base}.  If
  \code{recursive} is \code{TRUE}, the function searches \code{base} and
  all subdirectories.  A two-color IDAT files are pair of files with
  names ending in \code{_Red.idat} or \code{_Grn.idat}.

  If the \code{targets} argument is not \code{NULL} it is assumed it has
  a columned named \code{Basename}, and this is assumed to be pointing
  to the base name of a two color IDAT file, ie. a name that can be made
  into a real IDAT file by appending either \code{_Red.idat} or
  \code{_Grn.idat}.

  The type of methylation array is guess by looking at the number of
  probes in the IDAT files.
}
\value{
  An object of class \code{"RGChannelSet"} or
  \code{"RGChannelSetExtended"}.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\link{read.metharray}} for the workhorse function,
  \code{\link{read.metharray.sheet}} for reading a sample sheet and
  \code{\linkS4class{RGChannelSet}} for the output class.
}
\examples{
if(require(minfiData)) {

baseDir <- system.file("extdata", package = "minfiData")
RGset <- read.metharray.exp(file.path(baseDir, "5723646052"))

}
}
