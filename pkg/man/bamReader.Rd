\name{bamReader}
\alias{bamReader}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bamReader (function)
}
\description{
Static function which returns an instance of class bamReader.
}
\usage{
bamReader(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
Name of BAM file to be opened by returned bamReader object.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Creates an instance of class bamReader and returns it.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An instance of class bamReader. Before reading alignments, there must be a call to "bamOpen".
}
\references{
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
%%  ~~who you are~~
W. Kaisers
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
bamReader class.
}
\examples{
# bamReader Example:
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
bamOpen(r)
al<-GetNextAlignment(r)
bamClose(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bamReader}
