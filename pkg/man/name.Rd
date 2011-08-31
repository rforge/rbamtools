\name{name}
\alias{name}
\alias{name<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
name (bamAlignment)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Reading and writing accessor for read name.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Character string.
}
\references{
%% ~put references to the literature/web site here ~
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
%%  ~~who you are~~
Wolfgang Kaisers
}
\examples{
# name(bamAlignment) Example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
name(al)
name(al)<-"newName"
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Name}
