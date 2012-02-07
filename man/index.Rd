\name{index}
\alias{index}
\alias{index<-}
\title{index (bamAlignment)}
\description{Reading and writing Accessor for index attribute of bamAlignment objects.}
\usage{index(object)}
\arguments{
  \item{object}{An instance of class bamAlignment}
}
\value{Bool}
\author{Wolfgang Kaisers}
\examples{
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
idx_file <- system.file("extdata", "accepted_hits.bam.bai", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
align<-getNextAlign(reader)
index(reader)<-idx_file
index.initialized(reader)
coords<-as.integer(c(0,0,1874))
range<-bamRange(reader,coords)
align<-getNextAlign(range)
bamClose(reader)
}
\keyword{ index }
