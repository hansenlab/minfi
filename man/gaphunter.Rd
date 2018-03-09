\name{gaphunter}
\alias{gaphunter}
\title{
  Find gap signals in 450k data
}
\description{
  This function finds probes in the Illumina 450k Array for which calculated beta values
  cluster into distinct groups separated by a defined threshold. It identifies, for these
  \sQuote{gaps signals} the number of groups, the size of these groups, and the samples
  in each group.
}
\usage{
  gaphunter(object, threshold=0.05, keepOutliers=FALSE,
            outCutoff=0.01, verbose=TRUE)
}

\arguments{
  \item{object}{An object of class (Genomic)RatioSet, (Genomic)MethylSet, or matrix.
    If one of the first two, code{\link{getBeta}} is used to calculate beta values. If a
    matrix, must be one of beta values.}
  \item{threshold}{The difference in consecutive, ordered beta values
    that defines the presence of a gap signal. Defaults to 5 percent.}
  \item{keepOutliers}{Should outlier-driven gap signals be kept in the
    results? Defaults to \code{FALSE}}
  \item{outCutoff}{Value used to identify gap signals driven by outliers.
    Defined as the percentage of the total sample size; the sum of
    samples in all groups except the largest must exceed this number of
    samples in order for the probe to still be considered a gap signal.
    Defaults to 1 percent.}
 \item{verbose}{logical value. If \code{TRUE}, it writes some messages
   indicating progress. If \code{FALSE} nothing should be printed.}
}
\details{

  The function can calculate a beta matrix or utilize a user-supplied matrix of
  beta values.

  The function will idenfity probes with a gap in a beta signal greater than
  or equal to the defined \code{threshold}. These probes constitue an additional,
  dataset-specific subset of probes that merit special consideration due to their tendency
  to be driven by an underlying SNP or other genetic variant. In this manner, these probes
  can serve as surrogates for underlying genetic signal locally and/or in a broader
  (i.e. haplotype) context. Please see our upcoming manuscript for a detailed description
  of the utility of these probes.

  Outlier-driven gap signals are those in which the sum of the smaller group(s) does
  not exceed a certain percentage of the sample size, defined by the argument outCutoff.
}
\value{
  A list with three values,
  \item{proberesults}{A data frame listing, for each identified gap signal, the number
    of groups and the size of each group.}
  \item{sampleresults}{a matrix of dimemsions probes (rows) by samples (columns). Individuals
    are assigned numbers based onthe groups into which they cluster. Lower number groups
    indicate lower mean methylation values for the group. For example, individuals coded as
    \sQuote{1} will have a lower mean methylation value than those individuals coded as
    \sQuote{2}.}
  \item{algorithm}{A list detailing the arguments supplied to the function.}
}
\references{
  SV Andrews, C Ladd-Acosta, AP Feinberg, KD Hansen, MD Fallin.
  \emph{\sQuote{Gap hunting} to characterize clustered probe signals in Illumina methylation array data.}
  Epigenetics & Chromatin (2016) 9:56.
  doi:\href{http://dx.doi.org/10.1186/s13072-016-0107-z}{10.1186/s13072-016-0107-z}.
}
\author{
  Shan V. Andrews \email{sandre17@jhu.edu}.
}
\examples{
if(require(minfiData)) {
  gapres <- gaphunter(MsetEx.sub, threshold=0.3, keepOutliers=TRUE)
  #Note: the threshold argument is increased from the default value in this small example
  #dataset with 6 people to avoid the reporting of a large amount of probes as gap signals.
  #In a typical EWAS setting with hundreds of samples, the default arguments should be
  #sufficient.
}
}
