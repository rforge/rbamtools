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
    \item{\code{filename}:}{Object of class \code{"character"} ~~ }
    \item{\code{reader}:}{Object of class \code{"externalptr"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{bamClose}{\code{signature(object = "bamReader")}: ... }
    \item{bamOpen}{\code{signature(object = "bamReader")}: ... }
    \item{CreateIndex}{\code{signature(object = "bamReader")}: ... }
    \item{filename<-}{\code{signature(object = "bamReader")}: ... }
    \item{filename}{\code{signature(object = "bamReader")}: ... }
    \item{GetHeaderText}{\code{signature(object = "bamReader")}: ... }
    \item{GetNextAlignment}{\code{signature(object = "bamReader")}: ... }
    \item{GetReferenceCount}{\code{signature(object = "bamReader")}: ... }
    \item{GetReferenceData}{\code{signature(object = "bamReader")}: ... }
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
W.Kaisers
}
\note{
note
}
\seealso{
seealso
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("bamReader")
}
\keyword{classes}
