\name{QueryBases}
\alias{QueryBases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
QueryBases (bamAlignment)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Reading accessor for QueryBases attribute of bamAlignment.
}
\usage{
QueryBases(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An instance of bamAlignment class.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Returns input base-sequence of alignment.
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
Qualities
}
\examples{
# Qualities(bamAlignment)
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
QueryBases(al)
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{QueryBases}
