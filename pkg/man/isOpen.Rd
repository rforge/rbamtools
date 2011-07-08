\name{isOpen}
\alias{isOpen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
isOpen (bamReader, bamWriter)
}
\description{
Returns file-open status for bamReader and bamWriter objects.
}
\usage{
isOpen(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An instance of bamReader or bamWriter class
}
}
\details{
Opens BAM file.
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
\seealso{
bamOpen
}
\examples{
# isOpen Generic
# InsertSize (bamAlignment)
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
isOpen(r)
bamOpen(r)
isOpen(r)
# ...
bamClose(r)
isOpen(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{isOpen}
