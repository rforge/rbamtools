\name{bamSave}
\alias{bamSave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bamSave (bamWriter)
}
\description{
Saves bamAlignment to BAM file.
}
\usage{
bamSave(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An opened instance of class bamWriter
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
An instance of class bamAlignment
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
\examples{
# bamSave example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
w<-bamWriter("gapped_aligns_2.bam",GetReferenceData(r),samHeader(r))
bamOpen(w)
al<-GetNextAlignment(r)
bamSave(w,al)
bamClose(r)
bamClose(w)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bamSave }
