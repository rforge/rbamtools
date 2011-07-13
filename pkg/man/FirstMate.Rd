\name{FirstMate}
\alias{FirstMate}
\alias{FirstMate<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FirstMate (bamAlignment)
}
\description{
Reading and writing accessor function for (is)first-mate property of bamAlignment object.
}
\usage{
FirstMate(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An instance of class bamAlignment
}
}
\details{
Function can be used to read and write FirstMate status.
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
# FirstMate example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
FirstMate(al)
# TRUE or FALSE
FirstMate(al)<-TRUE
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{FirstMate}
