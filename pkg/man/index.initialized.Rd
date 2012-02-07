\name{index.initialized}
\alias{index.initialized}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{index.initialized(bamReader)}
\description{Checks for initialized index in bamReader.}
\usage{index.initialized(object)}
\arguments{
  \item{object}{An instance of class bamReader.}
}
\details{The function acutally checks, if the externalptr (index) slot has nil-value.}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx_file <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
index.initialized(reader)
index(reader)<-idx_file
index.initialized(reader)
bamClose(reader)
index.initialized(reader)
}
\keyword{index.initialized}
\keyword{bamOpen}
