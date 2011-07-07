\name{FailedQC}
\alias{FailedQC}
\alias{FailedQC<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FailedQC (bamAlignment) 
}
\description{
FailedQC denotes failed quality Control for read.
FailedQC method facilitates reading and writing access to FailedQC attribute of bamAlignemnt object.
}
\usage{
FailedQC(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An instance of class bamAlignment
}
}
\value{
bool
}
\references{
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009)
}
\author{
W.Kaisers
}
\examples{
# FailedQC(bamAlignment)
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
FailedQC(al)
FailedQC(al)<-TRUE
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{FailedQC}
