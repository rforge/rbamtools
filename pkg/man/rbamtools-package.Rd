\name{rbamtools-package}
\alias{rbamtools-package}
\alias{rbamtools}
\docType{package}
\title{
Reading and writing BAM files
}
\description{
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab rbamtools\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-06-24\cr
License: \tab GPL \cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
W.Kaisers
Maintainer: W.Kaisers <kaisers@med.uni-duesseldorf.de>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\keyword{ package }
\examples{
#############################
# rbamtools examples
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
#############################
# bamReader
r<-bamReader(bam_file)
filename(r)
bamOpen(r)
#############################
# bamAlignment
al<-GetNextAlignment(r)
name(al)
CigarData(al)
Qualities(al)
QueryBases(al)
AlignedBases(al)
#############################
# bamWriter
w<-bamWriter("test.bam",GetReferenceData(r),samHeader(r))
bamOpen(w)
bamSave(w,al)
bamClose(r)
bamClose(w)
#############################
}
