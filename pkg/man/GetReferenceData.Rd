\name{GetReferenceData}
\alias{GetReferenceData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GetReferenceData (bamReader)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns Reference data from BAM file.
}
\usage{
GetReferenceData(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An open instance of class bamReader.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Returns a data.frame with four columns (RefID, RefName, RefLength, RefHasAlignments). The function prints an error message and terminates when BAM files is not open.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
data.frame
}
\references{
%% ~put references to the literature/web site here ~
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College 2009
}
\author{
%%  ~~who you are~~
Wolfgang Kaisers
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
GetReferenceCount
}
\examples{
# GetReferenceData example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
GetReferenceData(r)
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GetReferenceData}
