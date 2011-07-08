\name{Duplicate}
\alias{Duplicate}
\alias{Duplicate<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Duplicate (bamAlignment)
}
\description{
Reading and writing Accessor for Duplicate attribute of bamAlignment objects.
}
\usage{
Duplicate(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An instance of class bamAlignment
}
}
\value{
Bool.
}
\references{
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
W.Kaisers
}
\examples{
# Duplicate
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
Duplicate(al)
Duplicate(al)<-TRUE
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Duplicate }
