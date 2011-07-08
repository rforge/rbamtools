\name{AlignedBases}
\alias{AlignedBases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AlignedBases (bamAlignment)
}
\description{
Returns Sequence of aligned bases.
}
\usage{
AlignedBases(object)
}
\arguments{
  \item{object}{
An instance of class bamAlignment
}
}
\details{
Returned Sequence contains 'N' characters for gapped aligns.
}
\value{
Character string
}
\references{
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
W.Kaisers
}
\seealso{
QueryBases
Qualities
}
\examples{
# BamAlignemnt.AlignedBases example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
AlignedBases(al)
bamClose(r)
# "CGATCTCCNNNNNNNNN...GCAGTTCTACC"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Aligned Bases}
