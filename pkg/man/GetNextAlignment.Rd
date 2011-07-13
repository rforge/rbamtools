\name{GetNextAlignment}
\alias{GetNextAlignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GetNextAlignment (bamReader)
}
\description{
Retrieves next Alignment from opened bamReader object.
}
\usage{
GetNextAlignment(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
A valid instance of class bamReader
}
}
\details{
Returns a new created object of type bamAlignment from bamReader object.
}
\seealso{
getNext
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
None.
}
\references{
%% ~put references to the literature/web site here ~
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
%%  ~~who you are~~
W.Kaisers
}
\examples{
# GetNextAlignment example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GetNextAlignment}
\keyword{bamReader}% __ONLY ONE__ keyword per line
