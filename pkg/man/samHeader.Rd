\name{samHeader}
\alias{samHeader}
\alias{samHeader<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
samHeader (bamReader, bamWriter)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Reading and writing access to SAM-header data in bamWriter and reading access in bamReader.
}
\usage{
samHeader(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
bamReader or bamWriter object.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Returns SAM-header data for bamReader and bamWriter. Can be used to replace SAM-header data in bamWriter.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
String.
}
\references{
%% ~put references to the literature/web site here ~
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
%%  ~~who you are~~
W.Kaisers
}
\examples{
# samHeader(bamReader) example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
samHeader(r)
bamClose(r)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{samHeader}
