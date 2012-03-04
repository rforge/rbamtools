\name{setVal}
\alias{setVal}
\title{setVal(object,members,values)}
\description{Setting values of data-members for object}
\usage{setVal(object,members,values)}
\arguments{
  \item{object}{An instance of class headerLine, refSeqDict, headerReadGroup or headerProgram }
  \item{members}{A (character) vector of member identifiers corresponding to data-members of the given object (e.g. PN for headerProgram)}
  \item{values}{Values that are written into data members}
}
\details{
The vectors members and values must have the same length.

Data members for class headerLine:
\tabular{ll}{
VN \tab Format version \cr
SO \tab Sorting order  \cr
}
Valid entries for SO are unknown (default), unsorted, queryname or coordinate.
Data members for class refSeqDict:
\tabular{ll}{
SN \tab Reference sequence name \cr
LN \tab Reference sequence length \cr
AS \tab Genome assembly identifier \cr
M5 \tab MD5 checksum of the sequence \cr
SP \tab Species \cr
UR \tab URI of the sequence \cr
}

Data members for class readGroup:
\tabular{ll}{
ID \tab Read Group identifier \cr
CN \tab Name of sequencing center \cr
DS \tab Description \cr
FO \tab Flow order \cr
KS \tab Nucleotides corresponding to key sequence of each read \cr
LB \tab Library \cr
PG \tab Programs used for processing the Read Group \cr
PI \tab Predicted median insert size \cr
PL \tab Sequencing Platform \cr
SM \tab Sample name \cr
}
Valid entries for PL are CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT or PACBIO

Data members for class headerProgram
\tabular{ll}{
ID \tab Program record identifier \cr
PN \tab Program name \cr
CL \tab Command line \cr
PP \tab Previous @PG-ID (Must match another ID in list) \cr
VN \tab Program Version \cr
}

}
\author{Wolfgang Kaisers}
\examples{
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
reader<-bamReader(bam_file)
header<-getHeader(reader)
hl<-getHeaderLine(header)
setVal(hl,"SO","coordinate")
bamClose(reader)
}
\keyword{setVal}
