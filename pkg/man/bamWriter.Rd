\name{bamWriter}
\alias{bamWriter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bamWriter (function)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns an instance of class bamWriter.
}
\usage{
bamWriter(filename, refData, samHeader)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
Name of BAM file to be opened
}
  \item{refData}{
%%     ~~Describe \code{refData} here~~
RefData (data.frame) which is written to BAM file.
}
  \item{samHeader}{
%%     ~~Describe \code{samHeader} here~~
samHeader (character string) which is written to BAM file.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Creates and returns a bamWriter object which has to be opened before writing.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An instance of class bamWriter.
}
\references{
%% ~put references to the literature/web site here ~
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
%%  ~~who you are~~
W. Kaisers
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
bamReader (function)
}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
w<-bamWriter("gapped_aligns_2.bam",GetReferenceData(r),samHeader(r))
bamOpen(w)
al<-GetNextAlignment(r)
bamSave(w,al)
bamClose(r)
bamClose(w)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bamWriter}
