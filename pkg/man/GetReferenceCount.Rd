\name{GetReferenceCount}
\alias{GetReferenceCount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GetReferenceCount (bamReader)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Retrieves number of References in opened BAM file. The function prints an error message and terminates when BAM file is not open. 
}
\usage{
GetReferenceCount(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An open instance of class bamReader (isOpen returns TRUE).
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Integer value.
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
GetReferenceData
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# GetReferenceCount example
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
GetReferenceCount(r)
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GetReferenceCount}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
