\name{read.450k.exp}
\alias{read.450k.exp}
\title{
  Reads an entire 450k experiment using a sample sheet
}
\description{
  Reads an entire 450k experiment using a sample sheet or (optionally)
  a target like data.frame.
}
\usage{
read.450k.exp(base, targets = NULL, extended = FALSE,
    recursive = FALSE, verbose = FALSE)
}
\arguments{
  \item{base}{The base directory.}
  \item{targets}{A targets \code{data.frame}, see details}
  \item{extended}{Should the output of the function be a
    \code{"RGChannelSetExtended"} (default is \code{"RGChannelSet"}).}
  \item{recursive}{Should the search be recursive (see details)}
  \item{verbose}{Should the function be verbose?}
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
}
\value{
  An object of class \code{"RGChannelSet"} or
  \code{"RGChannelSetExtended"}.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\link{read.450k}} for the workhorse function,
  \code{\link{read.450k.sheet}} for reading a sample sheet and
  \code{\linkS4class{RGChannelSet}} for the output class.
}
\examples{
if(require(minfiData)) {

baseDir <- system.file("extdata", package = "minfiData")
RGset <- read.450k.exp(file.path(baseDir, "5723646052"))

}
}
