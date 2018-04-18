\name{subsetByLoci}
\alias{subsetByLoci}
\title{
  Subset an RGChannelset by CpG loci.
}
\description{
  Subset an RGChannelSet by CpG loci.
}
\usage{
subsetByLoci(rgSet, includeLoci = NULL, excludeLoci = NULL,
             keepControls = TRUE, keepSnps = TRUE)
}
\arguments{
  \item{rgSet}{An object of class \code{RGChannelSet} (or
    \code{RGChannelSetExtended}).}
  \item{includeLoci}{A character vector of CpG identifiers which should
    be kept.}
  \item{excludeLoci}{A character vector of CpG identifiers which should
    be excluded.}
  \item{keepControls}{Should control probes be kept?}
  \item{keepSnps}{Should SNP probes be kept?}
}
\details{
  This task is non-trivial because an \code{RGChannelSet} is indexed by
  probe position on the array, not by loci name.
}
\value{
  An object of class \code{RGChannelSet}, which some probes removed.
}
\examples{
if(require(minfiData)) {
  loci <- c("cg00050873", "cg00212031", "cg00213748", "cg00214611")
  subsetByLoci(RGsetEx.sub, includeLoci = loci)
  subsetByLoci(RGsetEx.sub, excludeLoci = loci)
}
}
