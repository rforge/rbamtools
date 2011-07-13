\name{CreateIndex}
\alias{CreateIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CreateIndex (bamReader)
}
\description{
Creates index file ("filename.bam.bai") and returns TRUE on success.
}
\usage{
CreateIndex(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An open instance of bamReader class. When bamReader is not open, an error occurs.
}
}
\value{
Bool.
}
\references{
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
Wolfgang Kaisers
}
\examples{
# CreateIndex example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
CreateIndex(r)
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{CreateIndex}
