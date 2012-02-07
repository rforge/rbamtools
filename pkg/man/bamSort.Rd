\name{bamSort}
\alias{bamSort}
\title{bamSort(bamAlign)}
\description{Sorts a given BAM-file into a second BAM-file.}
\usage{bamSort(object,prefix,byName,maxmem)}
\arguments{
  \item{object}{An Instance of class bamReader}
  \item{prefix}{Character String}
  \item{byName}{Logical}
  \item{maxmem}{Approximate maximum memory usage}
}
\details{\verb{BAM-files must be sorted before an index can be created (Otherwise an error it thrown). During this routine, some content of the source file is copied into the RAM, sorted and written into a temporary file. The number and size of temporary files depend on the size of the source file and the utilized RAM. Finally, the sorted files are merged into a single target file. The name of the target file is [prefix].bam. When byName is FALSE (standard), the BAM-file is sorted by refID and position. When byName is TRUE, the BAM-file is sorted by name.}}
\value{None}
\author{Wolfgang Kaisers}

\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
isOpen(reader)
bamSort(reader,"sorted")
bamClose(reader)
}
\keyword{bamReader}
\keyword{bamSort}
