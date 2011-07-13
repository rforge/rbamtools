\name{Qualities}
\alias{Qualities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Qualities (bamAlignment)
}
\description{
Returns string with ASCII-coded Quality scores.
}
\usage{
Qualities(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An instance of class bamAlignment.
}
}
\details{
Besides the nucleotide sequence of a read, an alignment program usually reports the quality scores for each nucleotide. These quality scores are ASCII-coded into a character sequence. The quality scores should and the nucleotide sequence (QueryBases) should have the same length.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A character string.
}
\references{
%% ~put references to the literature/web site here ~
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
%%  ~~who you are~~
W.Kaisers
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
AlignedBases
QueryBases
}
\examples{
# Qualities(bamAlignment)
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
Qualities(al)
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Qualities}
