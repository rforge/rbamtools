\name{isValid}
\alias{isValid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
isValid (bamAlignment)
}
\description{
Returns validity status of bamAlignment object.
}
\usage{
isValid(object)
}
\arguments{
  \item{object}{
An instance of class bamAlignment
}
}
\details{
A successful call of GetNextAlignment returns a valid bamAlignment object. Calling GetNextAlignment when read position in BAM file is past the last record retrieves an invalid bamAlignment object.
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
\examples{
# BamAlignemnt.isValid example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
while(isValid(al))
{
	# ...
	al<-GetNextAlignment(r)
}
# ...
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{isValid}
