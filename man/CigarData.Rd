\name{CigarData}
\alias{CigarData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CigarData (bamAlignment)
}
\description{
Reading accessor for cigar-data.
}
\usage{
CigarData(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An instance of class bamAlignment
}
}
\value{
data.frame with four columns: RefID, RefName, RefLength, RefHasAlignments
}
\references{
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009)
}
\author{
W.Kaisers
}
\examples{
#CigarData example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
cigar<-CigarData(al)
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Cigar data}
