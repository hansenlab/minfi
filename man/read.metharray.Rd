\name{read.metharray}
\alias{read.metharray}
\title{
  Parsing IDAT files from Illumina methylation arrays.
}
\description{
  Parsing IDAT files from Illumina methylation arrays.
}
\usage{
read.metharray(basenames, extended = FALSE, verbose = FALSE, force = FALSE)
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
  \item{force}{Should reading different size IDAT files be forced? See
    Details.}
}
\details{
  The type of methylation array is guess by looking at the number of
  probes in the IDAT files.

  We have seen IDAT files from the same array, but with different number
  of probes in the wild.  Specifically this is the case for early access
  EPIC arrays which have fewer probes than final release EPIC arrays.
  It is possible to combine IDAT files from the same inferred array, but
  with different number of probes, into the same \code{RGChannelSet} by
  setting \code{force=TRUE}.  The output object will have the same
  number of probes as the smallest array being parsed; effectively
  removing probes which could have been analyzed.
}
\value{
  An object of class \code{RGChannelSet} or
  \code{RGChannelSetExtended}. 
}
\author{
  Kasper Daniel Hansen\email{khansen@jhsph.edu}.
}
\seealso{
  \code{\link{read.metharray.exp}} for a convenience function for
  reading an experiment, \code{\link{read.metharray.sheet}} for
  reading a sample sheet and \code{\linkS4class{RGChannelSet}} for the
  output class. 
}
\examples{
if(require(minfiData)) {

baseDir <- system.file("extdata", package = "minfiData")
RGSet <- read.metharray(file.path(baseDir, "5723646052",
    c("5723646052_R02C02", "5723646052_R04C01")))

}
}
