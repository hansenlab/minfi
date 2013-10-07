\name{blockFinder}
\alias{blockFinder}
\title{
  Finding blocks of methylation differences for Illumina methylation arrays
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
blockFinder(object, design, coef = 2, what = c("Beta", "M"),
            cluster = NULL, cutoff = NULL, pickCutoff = FALSE,
            smooth = TRUE, smoothFunction = locfitByCluster,
            B = 1000, verbose = TRUE, bpSpan = 2.5 * 10^5, ...)
}
\arguments{
  \item{object}{}
  \item{design}{}
  \item{coef}{}
  \item{what}{}
  \item{cluster}{}
  \item{cutoff}{}
  \item{pickCutoff}{}
  \item{smooth}{}
  \item{smoothFunction}{}
  \item{B}{}
  \item{verbose}{}
  \item{bpSpan}{}
  \item{\dots}{}
}
\details{
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
}
\note{
}
\seealso{
  \code{\link{bumphunter}}
}
\examples{
if(require(minfiData)) {

}
}
