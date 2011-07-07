\name{bamAlignment-class}
\Rdversion{1.1}
\docType{class}
\alias{bamAlignment-class}
\alias{AlignedBases,bamAlignment-method}
\alias{CigarData,bamAlignment-method}
\alias{Duplicate<-,bamAlignment-method}
\alias{Duplicate,bamAlignment-method}
\alias{FailedQC<-,bamAlignment-method}
\alias{FailedQC,bamAlignment-method}
\alias{FirstMate<-,bamAlignment-method}
\alias{FirstMate,bamAlignment-method}
\alias{getNext,bamAlignment-method}
\alias{initialize,bamAlignment-method}
\alias{InsertSize,bamAlignment-method}
\alias{isValid,bamAlignment-method}
\alias{Mapped<-,bamAlignment-method}
\alias{Mapped,bamAlignment-method}
\alias{MapQuality,bamAlignment-method}
\alias{MateMapped<-,bamAlignment-method}
\alias{MateMapped,bamAlignment-method}
\alias{MatePosition,bamAlignment-method}
\alias{MateRefID,bamAlignment-method}
\alias{MateReverseStrand<-,bamAlignment-method}
\alias{MateReverseStrand,bamAlignment-method}
\alias{name,bamAlignment-method}
\alias{name<-,bamAlignment-method}
\alias{Paired<-,bamAlignment-method}
\alias{Paired,bamAlignment-method}
\alias{position,bamAlignment-method}
\alias{PrimaryAlignment<-,bamAlignment-method}
\alias{PrimaryAlignment,bamAlignment-method}
\alias{ProperPair<-,bamAlignment-method}
\alias{ProperPair,bamAlignment-method}
\alias{Qualities,bamAlignment-method}
\alias{QueryBases,bamAlignment-method}
\alias{refID,bamAlignment-method}
\alias{ReverseStrand<-,bamAlignment-method}
\alias{ReverseStrand,bamAlignment-method}
\alias{SecondMate<-,bamAlignment-method}
\alias{SecondMate,bamAlignment-method}

\title{S4 Class "bamAlignment"}
\description{
Container for all Alignment information usually read from BAM-file. Can also be saved to open BAM-file with bamWriter.
}
\section{Objects from the Class}{
Usually retrieved with GetNextAlignment on open bamReader Object. 
}
\section{Slots}{
  \describe{
    \item{\code{bamAlignment}:}{Object of class \code{"externalptr"} wich points to external BamAlignment object. }
    \item{\code{valid}:}{Boolean indicating that Object contains valid information. Invalid bamAlignment is returned when GetNextAlignment is called after the last Alignment set has been retrieved.}
    \item{\code{reader}:}{Object of class \code{"externalptr"} which points to external bamReader object.}
}}
\section{Methods}{
  \describe{
    \item{AlignedBases}{\code{signature(object = "bamAlignment")}: ... }
    \item{CigarData}{\code{signature(object = "bamAlignment")}: ... }
    \item{Duplicate<-}{\code{signature(object = "bamAlignment")}: ... }
    \item{Duplicate}{\code{signature(object = "bamAlignment")}: ... }
    \item{FailedQC<-}{\code{signature(object = "bamAlignment")}: ... }
    \item{FailedQC}{\code{signature(object = "bamAlignment")}: ... }
    \item{FirstMate<-}{\code{signature(object = "bamAlignment")}: ... }
    \item{FirstMate}{\code{signature(object = "bamAlignment")}: ... }
    \item{getNext}{\code{signature(object = "bamAlignment")}: ... }
    \item{initialize}{\code{signature(.Object = "bamAlignment")}: ... }
    \item{InsertSize}{\code{signature(object = "bamAlignment")}: ... }
    \item{isValid}{\code{signature(object="bamAlignment")}: ... }
    \item{Mapped<-}{\code{signature(object = "bamAlignment")}: ... }
    \item{Mapped}{\code{signature(object = "bamAlignment")}: ... }
    \item{MapQuality}{\code{signature(object = "bamAlignment")}: ... }
    \item{MateMapped<-}{\code{signature(object = "bamAlignment")}: ... }
    \item{MateMapped}{\code{signature(object = "bamAlignment")}: ... }
    \item{MatePosition}{\code{signature(object = "bamAlignment")}: ... }
    \item{MateRefID}{\code{signature(object = "bamAlignment")}: ... }
    \item{MateReverseStrand<-}{\code{signature(object = "bamAlignment")}: ... }
    \item{MateReverseStrand}{\code{signature(object = "bamAlignment")}: ... }
    \item{name}{\code{signature(object = "bamAlignment")}: ... }
    \item{Paired<-}{\code{signature(object = "bamAlignment")}: ... }
    \item{Paired}{\code{signature(object = "bamAlignment")}: ... }
    \item{position}{\code{signature(object = "bamAlignment")}: ... }
    \item{PrimaryAlignment<-}{\code{signature(object = "bamAlignment")}: ... }
    \item{PrimaryAlignment}{\code{signature(object = "bamAlignment")}: ... }
    \item{ProperPair<-}{\code{signature(object = "bamAlignment")}: ... }
    \item{ProperPair}{\code{signature(object = "bamAlignment")}: ... }
    \item{Qualities}{\code{signature(object = "bamAlignment")}: ... }
    \item{QueryBases}{\code{signature(object = "bamAlignment")}: ... }
    \item{refID}{\code{signature(object = "bamAlignment")}: ... }
    \item{ReverseStrand<-}{\code{signature(object = "bamAlignment")}: ... }
    \item{ReverseStrand}{\code{signature(object = "bamAlignment")}: ... }
    \item{SecondMate<-}{\code{signature(object = "bamAlignment")}: ... }
    \item{SecondMate}{\code{signature(object = "bamAlignment")}: ... }
	 }
}
\references{
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009)
}
\author{
W.Kaisers
}
\examples{
# bamAlignment-class examples
library(rbamtools)
bam_file <- system.file("extdata", "accepted_hits.bam", package="rbamtools")
r<-bamReader(bam_file)
filename(r)
bamOpen(r)
al<-GetNextAlignment(r)
bamClose(r)
}
\keyword{classes}
