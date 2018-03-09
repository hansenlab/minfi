\name{preprocessRaw}
\alias{preprocessRaw}
\title{
  Creation of a MethylSet without normalization
}
\description{
  Converts the Red/Green channel for an Illumina methylation array into
  methylation signal, without using any normalization.
}
\usage{
preprocessRaw(rgSet)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet}.}
}
\details{
  This function takes the Red and the Green channel of an Illumina
  methylation array, together with its associated manifest object and
  converts it into a \code{MethylSet} containing the methylated and
  unmethylated signal.
}
\value{
  An object of class \code{MethylSet}
}
\author{
  Kasper Daniel Hansen\email{khansen@jhsph.edu}.
}
\seealso{
  \code{\linkS4class{RGChannelSet}} and \code{\linkS4class{MethylSet}}
  as well as \code{\linkS4class{IlluminaMethylationManifest}}.
}
\examples{
if (require(minfiData)) {

dat <- preprocessRaw(RGsetEx)
slot(name="preprocessMethod", dat)[1]

}
}
