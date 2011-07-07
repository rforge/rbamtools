\name{bamClose}
\alias{bamClose}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bamClose: Generic function for closing file connection for bamReader and bamWriter
}
\description{
Closes BAM file for bamReader and bamWriter objects.
}
\usage{
bamClose(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An instance of class bamReader or bamWriter.
}
}
\references{ Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009)
}
\author{
W.Kaisers
}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
filename(r)
bamOpen(r)
al<-GetNextAlignment(r)
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bamClose}
\keyword{bamOpen}
