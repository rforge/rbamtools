\name{bamOpen}
\alias{bamOpen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bamOpen (bamReader, bamWriter)
}
\description{
Generic function for bamReader and bamWriter for opening BAM file.
}
\usage{
bamOpen(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An instance of bamReader or bamWriter class
}
  \item{...}{
%%     ~~Describe \code{\dots} here~~
Optional name of BAM Index file for bamReader.
}
}
\value{
None.
}
\references{
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
Wolfgang Kaisers
}
\seealso{
bamClose
}
\examples{
# bamOpen example
# BamAlignemnt.AlignedBases example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
# ...
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bamOpen }
