\name{Paired}
\alias{Paired}
\alias{Paired<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Paired (bamAlignment)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns information whether alignment is part of paired-end read.
}
\usage{
Paired(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An instance of class bamAlignment.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Paired
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Bool.
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
MateMapped
MateRefID
MatePosition
MateReverseStrand
}
\examples{
# Paired(bamAlignment)
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
Paired(al)
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Paired}
