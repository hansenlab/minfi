\name{cpgCollapse}
\alias{cpgCollapse}
\title{
  Collapse methylation annotation for block finding
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cpgCollapse(object, what = c("Beta", "M"), maxGap = 500,
            blockMaxGap = 2.5 * 10^5, maxClusterWidth = 1500,
            dataSummary = colMeans, na.rm = FALSE,
            returnBlockInfo = TRUE, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{}
  \item{what}{}
  \item{maxGap}{}
  \item{blockMaxGap}{}
  \item{maxClusterWidth}{}
  \item{dataSummary}{}
  \item{na.rm}{}
  \item{returnBlockInfo}{}
  \item{verbose}{}
  \item{\dots}{}
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
}
\author{
}
\seealso{
  \code{\link{blockFinder}}
}
\examples{
if(require(minfiData)) {

}
}
