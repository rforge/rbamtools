\name{getNext}
\alias{getNext}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getNext (bamAlignment)
}
\description{
Fills next alignment from BAM file into bamAlignment object without creating a new instance.
}
\usage{
getNext(.Object)
}
\arguments{
  \item{.Object}{
An instance of class bamAlignment
}
}
\value{
None.
}
\author{
W.Kaisers
}
\seealso{
GetNextAlignment
}
\examples{
# BamAlignemnt.AlignedBases example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
name(al)
getNext(al)
name(al)
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{getNext}
