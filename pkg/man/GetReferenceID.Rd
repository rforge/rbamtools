\name{GetReferenceID}
\alias{GetReferenceID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GetReferenceID (bamReader)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
GetReferenceID returns refID from BAM file for given RefName.
}
\usage{
GetReferenceID(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An open instance of class bamReader.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
The name of the reference to be searched for.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
GetReferenceID returns the refID for a given RefName in a BAM file (represented by bamReader). The function prints an error message and terminates when bamReader is not open.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Integer.
}
\references{
%% ~put references to the literature/web site here ~
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
%%  ~~who you are~~
Wolfgang Kaisers
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
GetReferenceData, GetReferenceCount
}
\examples{
# GetReferenceID example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
GetReferenceID(r,"NT_004350")
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GetReferenceID}
