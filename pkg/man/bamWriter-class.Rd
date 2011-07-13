\name{bamWriter-class}
\Rdversion{1.1}
\docType{class}
\alias{bamWriter-class}
\alias{initialize-methods}
\alias{filename<-,bamWriter-method}
\alias{filename,bamWriter-method}
\alias{initialize,bamWriter-method}
\alias{refData<-,bamWriter-method}
\alias{refData,bamWriter-method}
\alias{samHeader<-,bamWriter-method}
\alias{samHeader,bamWriter-method}

\title{Class "bamWriter" Class for writing bamAlignments to BAM files}
\description{
bamWriter: saves bamAlignemnts to BAM files
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bamWriter", filename, refData, samHeader)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{filename}:}{Object of class \code{"character"} ~~ }
    \item{\code{writer}:}{Object of class \code{"externalptr"} ~~ }
    \item{\code{ref_data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{sam_header}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{bamClose}{\code{signature(object = "bamWriter")}: ... }
    \item{bamOpen}{\code{signature(object = "bamWriter")}: ... }
    \item{bamSave}{\code{signature(object = "bamWriter")}: ... }
    \item{filename<-}{\code{signature(object = "bamWriter")}: ... }
    \item{filename}{\code{signature(object = "bamWriter")}: ... }
    \item{initialize}{\code{signature(.Object = "bamWriter")}: ... }
    \item{refData<-}{\code{signature(object = "bamWriter")}: ... }
    \item{refData}{\code{signature(object = "bamWriter")}: ... }
    \item{samHeader<-}{\code{signature(object = "bamWriter")}: ... }
    \item{samHeader}{\code{signature(object = "bamWriter")}: ... }
	 }
}
\references{
Derek Barnett, Michael Stromberg, Marth Lab, Department of Biology, Boston College (2009).
}
\author{
Wolfgang Kaisers
}
\examples{
showClass("bamWriter")
}
\keyword{classes}
