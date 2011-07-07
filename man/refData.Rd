\name{refData}
\alias{refData}
\alias{refData<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
refData (bamWriter)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Reads and writes reference data in bamWriter.
}
\usage{
refData(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The data.frame passed from and to refData has four columns: RefID, RefName, RefLength, RefHasAlignments.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data.frame with four columns.
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
refData
}
\examples{
# refData example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
w<-bamWriter("corrected.bam",GetReferenceData(r),samHeader(r))
ref<-refData(w)
# ...
refData(w)<-ref
bamOpen(w)
# ...
bamClose(w)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{refData}
