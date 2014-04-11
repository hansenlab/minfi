\name{read.450k}
\alias{read.450k}
\title{
  Parsing IDAT files from Illumina methylation arrays.
}
\description{
  Parsing IDAT files from Illumina methylation arrays.
}
\usage{
read.450k(basenames, extended = FALSE, verbose = FALSE)
}
\arguments{
  \item{basenames}{The basenames or filenames of the IDAT files.  By
    basenames we mean the filename without the ending \code{_Grn.idat}
    or \code{_Red.idat} (such that each sample occur once).  By
    filenames we mean filenames including \code{_Grn.idat}
    or \code{_Red.idat} (but only one of the colors)}
  \item{extended}{Should a \code{RGChannelSet} or a
    \code{RGChannelSetExtended} be returned. }
  \item{verbose}{Should the function be verbose?}
}
\value{
  An object of class \code{RGChannelSet} or
  \code{RGChannelSetExtended}. 
}
\author{
  Kasper Daniel Hansen\email{khansen@jhsph.edu}.
}
\seealso{
  \code{\link{read.450k.exp}} for a convenience function for
  reading an experiment, \code{\link{read.450k.sheet}} for
  reading a sample sheet and \code{\linkS4class{RGChannelSet}} for the
  output class. 
}
\examples{
if(require(minfiData)) {

baseDir <- system.file("extdata", package = "minfiData")
RGSet <- read.450k(file.path(baseDir, "5723646052",
    c("5723646052_R02C02", "5723646052_R04C01")))

}
}
