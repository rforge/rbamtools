\name{Rewind}
\alias{Rewind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rewind (bamReader)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Resets the reading position of bamReader to beginning of file. The next call to GetNextAlignment then retrieves the first alignment in the opened BAM-file.
}
\usage{
Rewind(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
None.
}
\references{
%% ~put references to the literature/web site here ~
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
%%  ~~who you are~~
W.Kaisers
}
\note{
%%  ~~further notes~~
Rewind
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
Rewind
}
\examples{
# Rewind(bamReader)
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
position(GetNextAlignment(r))
Rewind(r)
position(GetNextAlignment(r))
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Rewind}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
