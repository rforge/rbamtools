\name{rbamtools-package}
\alias{rbamtools-package}
\alias{rbamtools}
\docType{package}
\title{
Reading and writing BAM files
}
\description{
The package contains functions for accessing contents of BAM (BGZF compressed SAM=Sequence Alignment/Map format) files.
}
\details{
\tabular{ll}{
Package: \tab rbamtools\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-06-20\cr
License: \tab GPL \cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}
Opening BAM files with bamReader objects. Sequential reading file content into bamAlignment objects. Accessing alignment content via accessor functions. Writing BAM files with bamWriter objects.
}
\author{
W.Kaisers
Maintainer: W.Kaisers <kaisers@med.uni-duesseldorf.de>
}
\references{
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
