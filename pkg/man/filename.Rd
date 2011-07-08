\name{filename}
\alias{filename}
\alias{filename<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
filename (bamReader, bamWriter)
}
\description{
Reading and writing accessor for filename attribute of bamReader and bamWriter objects.
Filename identifies the BAM file to wich bamReader and bamWriter are connected.
}
\usage{
filename(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An instance of class bamReader or bamWriter.
}
}

\value{
String containing current filename of bamReader object (for reading attribute access).
}
\references{
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
W.Kaisers
}
\examples{
# filename example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
filename(r)
filename(r)<-"new_file_name.bam"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{filename}
