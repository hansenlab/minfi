\name{bumphunter-methods}
\docType{methods}
\alias{bumphunter,GenomicRatioSet-method}

\title{Methods for function \code{bumphunter} in Package \pkg{minfi}}

\description{
  Estimate regions for which a genomic profile deviates from its baseline
  value. Originally implemented to detect differentially methylated
  genomic regions between two populations, but can be applied to any CpG-level 
  coefficient of interest.
  
}
\usage{
\S4method{bumphunter}{GenomicRatioSet}(object, design, cluster=NULL,
          coef=2, cutoff=NULL, cutoffQ = 0.99,
          maxGap=500, smooth=FALSE, smoothFunction=loessByCluster,
          useWeights=FALSE, B=1000, verbose=TRUE, 
          type = c("M", "Beta"), \dots)
}

\arguments{
  \item{object}{An object of class GenomicRatioSet.}
  \item{design}{Design matrix with rows representing samples and columns
    representing covariates. Regression is applied to each row of mat.}
  \item{cluster}{The clusters of locations that are to be analyzed
    together.  In the case of microarrays, the clusters are many times
    supplied by the manufacturer. If not available the function
    \code{\link{clusterMaker}} can be used to cluster nearby locations.}
  \item{coef}{An integer denoting the column of the design matrix
    containing the covariate of interest. The hunt for bumps will be
    only be done for the estimate of this coefficient.}
  \item{cutoff}{A numeric value. Values of the estimate of the genomic
    profile above the cutoff or below the negative of the cutoff will be
    used as candidate regions.  It is possible to give two separate
    values (upper and lower bounds).  If one value is given, the lower
    bound is minus the value.}
  \item{cutoffQ}{The quantile used for picking the cutoff using the
    permutation distribution.}
  \item{maxGap}{ If cluster is not provided this maximum location gap
    will be used to define cluster via the  \code{\link{clusterMaker}}
  function.}
  \item{smooth}{
    A logical value. If TRUE the estimated profile will be smoothed with the
    smoother defined by \code{smoothFunction}}
  \item{smoothFunction}{
    A function to be used for smoothing the estimate of the genomic
    profile. Two functions are provided by the package: \code{loessByCluster}
    and \code{runmedByCluster}.}
  \item{useWeights}{
    A logical value. If \code{TRUE} then the standard errors of the
   point-wise estimates of the profile function will be used as weights
   in the loess smoother \code{loessByCluster}. If the 
   \code{runmedByCluster} smoother is used this argument is ignored.}
  \item{B}{An integer denoting the number of resamples to use when computing
    null distributions: 1000 is recommended but can be quite slow if
    this process is not parallelized.}
  \item{verbose}{
    logical value. If \code{TRUE}, it writes out some messages
    indicating progress. If \code{FALSE} nothing should be printed.
  }
  \item{type}{
    Should bumphunting be performed on M-values ("M") or Beta values
    ("Beta")?}
  \item{\dots}{
    further arguments to be passed to the smoother functions.
}
}
\details{

This function performs the bumphunting approach described by Jaffe et
al. International Journal of Epidemiology (2012). The main output is a
table of candidate regions with permutation-based family-wide error
rates (FWER) and p-values assigned.

The general idea is that for each genomic location we have a value for
several individuals. We also have covariates for each individual and
perform regression. This gives us one estimate of the coefficient of
interest (a common example is case versus control). These estimates are
then (optionally) smoothed. The smoothing occurs in clusters of
locations that are \sQuote{close enough}. This gives us an estimate of a
genomic profile that is 0 when uninteresting. We then take values above
(in absolute value) \code{cutoff} as candidate regions. Permutations are
then performed to create null distributions for the candidate regions.

Uncertainty is assessed via permutations. Each of the B permutations
will produce an estimated \sQuote{null profile} from which we can define
\sQuote{null candidate regions}. For each observed candidate region we
determine how many null regions are \sQuote{more extreme} (longer and
higher average value). The \sQuote{p.value} is the percent of candidate
regions obtained from the permutations that as extreme as the observed
region. These p-values should be interpreted with care as the
theoretical proporties are not well understood.  The \sQuote{fwer} is
the proportion of permutations that had at least on regio as extreme as
the observed region. We compute p.values and FWER for the area of the
regions (as opposed to length and value as a pair) as well.

Parallelization is implemented through the foreach package.
}

\value{
An object of class \code{bumps} with the following components:
\item{tab}{The table with candidate regions and annotation for these.}
\item{coef}{The single loci coefficients.}
\item{fitted}{The estimated genomic profile used to determine the
  regions.}
\item{pvaluesMarginal}{marginal p-value for each genomic location.}
\item{null}{The null distribution.}
\item{algorithm}{details on the algorithm.}
}

\references{
Jaffe AE, Murakami P, Lee H, Leek JT, Fallin MD, Feinberg AP, Irizarry
RA (2012) Bump hunting to identify differentially methylated regions in
epigenetic epidemiology studies. International Journal of Epidemiology
41(1):200-9.
}

\author{
  Rafael A. Irizarry, Martin J. Aryee and Kasper D. Hansen
}


\examples{
if(require(minfiData)) {
  gmSet <- preprocessQuantile(MsetEx)
  design <- model.matrix(~ gmSet$status)
  bumps <- bumphunter(gmSet, design = design, B = 0,
                      type = "Beta", cutoff = 0.25)
}
}


