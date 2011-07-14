\name{bamReader-class}
\Rdversion{1.1}
\docType{class}
\alias{bamReader-class}
\alias{CreateIndex,bamReader-method}
\alias{filename<-,bamReader-method}
\alias{filename,bamReader-method}
\alias{samHeader,bamReader-method}
\alias{GetNextAlignment,bamReader-method}
\alias{GetReferenceCount,bamReader-method}
\alias{GetReferenceData,bamReader-method}
\alias{samHeader,bamReader-method}
\alias{GetReferenceID,bamReader-method}
\alias{initialize,bamReader-method}
\alias{Rewind,bamReader-method}

\title{Class "bamReader"}
\description{
bamReader provides Methods for opening and navigating BAM-files and to retrieve bamAlignment objects
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bamReader", filename)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{filename}:}{Object of class \code{"character"}:\tabName of the (connected) BAM-file. }
    \item{\code{reader}:}{Object of class \code{"externalptr"}:\tabPointer to external C++ Object of class BamReader }
  }
}
\section{Methods}{
  \describe{
    \item{bamClose}{\code{signature(object = "bamReader")}: Close opened BAM file }
    \item{bamOpen}{\code{signature(object = "bamReader")}: Open BAM file }
    \item{CreateIndex}{\code{signature(object = "bamReader")}: Create index file (.bam.bai) file for BAM-file. }
    \item{filename<-}{\code{signature(object = "bamReader")}: Setting name of BAM-file }
    \item{filename}{\code{signature(object = "bamReader")}: Reading name of BAM-file }
    \item{GetHeaderText}{\code{signature(object = "bamReader")}: Get SAM-Header (as text) from BAM-file }
    \item{GetNextAlignment}{\code{signature(object = "bamReader")}: Reading next alignment from BAM-file }
    \item{GetReferenceCount}{\code{signature(object = "bamReader")}: Reading number of reference sequences in alignment }
    \item{GetReferenceData}{\code{signature(object = "bamReader")}: Reading reference data as data.frame }
    \item{samHeader}{\code{signature(object = "bamReader")}: ... }
    \item{GetReferenceID}{\code{signature(object = "bamReader")}: ... }
    \item{initialize}{\code{signature(.Object = "bamReader")}: ... }
    \item{Rewind}{\code{signature(object = "bamReader")}: ... }
	 }
}
\references{
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
Wolfgang Kaisers
}
\examples{
showClass("bamReader")
}
\keyword{classes}
