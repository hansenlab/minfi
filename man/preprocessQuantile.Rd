\name{preprocessQuantile}
\alias{preprocessQuantile}
\title{
  Stratified quantile normalization for an Illumina methylation array.
}
\description{
  Stratified quantile normalization for Illumina amethylation arrays.
}
\usage{
preprocessQuantile(object, fixOutliers = TRUE, removeBadSamples = FALSE,
                   badSampleCutoff = 10.5, quantileNormalize = TRUE,
                   stratified = TRUE, mergeManifest = FALSE, sex = NULL,
                   verbose = TRUE)
}
\arguments{
  \item{object}{An object of class \code{RGChannelSet} or
    \code{[Genomic]MethylSet}.}
  \item{fixOutliers}{}
  \item{removeBadSamples}{}
  \item{badSampleCutoff}{}
  \item{quantileNormalize}{}
  \item{stratified}{}
  \item{mergeManifest}{}
  \item{sex}{}
  \item{verbose}{Should the function be verbose?}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Rafael A. Irizarry
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
if(require(minfiData)) {

}
}
