# #####################################################################################################################
# Generics definitions for bamReader and bamWriter
setGeneric("filename", function(object) standardGeneric("filename"))
setGeneric("filename<-", function(object,value) standardGeneric("filename<-"))
setGeneric("bamOpen", function(object,...) standardGeneric("bamOpen"))
setGeneric("bamClose", function(object) standardGeneric("bamClose"))
setGeneric("isOpen",function(object) standardGeneric("isOpen"))
setGeneric("samHeader", function(object) standardGeneric("samHeader"))
setGeneric("samHeader<-", function(object,value) standardGeneric("samHeader<-"))

# #####################################################################################################################
# bamReader
bamReader<-function(filename) { return(new("bamReader",filename=filename))}

setClass("bamReader", representation(filename="character",reader="externalptr"),
			validity=function(object){
				return(TRUE)
			}
)

setMethod(f="initialize", signature="bamReader",
		definition=function(.Object,filename){
			.Object@filename<-filename
			.Object@reader=.Call("get_bam_reader")
			return(.Object)
		}
)


# generic in header
setMethod("filename", "bamReader", function(object) return(object@filename))
# generic in header
setReplaceMethod(f="filename", signature="bamReader",
   definition=function(object,value){
       if(!is.character(value))
           stop("class bamReader, filename setter: value must be character")
       object@filename<-value
       return(object)
   }
)

# generic in header
setMethod(f="bamOpen",signature="bamReader", definition=function(object,index_file="")
				{
					if(!file.exists(object@filename))
					{
						stop(paste("bamReader.bamOpen: file",object@filename,"does not exist!"))
					}
					.Call("bam_reader_open",object@reader,object@filename,index_file)
					return(invisible())
				})

# generic in header
setMethod(f="bamClose",signature="bamReader",definition=function(object) {.Call("bam_reader_close",object@reader);return(invisible())})
# generic in header
setMethod("isOpen",signature="bamReader",definition=function(object) {.Call("bam_reader_isOpen",object@reader)})

#setGeneric("GetHeaderText", function(object) standardGeneric("GetHeaderText"))
# generic in header
setMethod(f="samHeader",signature="bamReader",definition=function(object) {.Call("bam_reader_GetHeaderText",object@reader)})
setGeneric("GetReferenceCount",function(object) standardGeneric("GetReferenceCount"))
setMethod(f="GetReferenceCount",signature="bamReader",definition=function(object) {.Call("bam_reader_GetReferenceCount",object@reader)})
setGeneric("GetReferenceData",function(object) standardGeneric("GetReferenceData"))
setMethod(f="GetReferenceData",signature="bamReader",definition=function(object) {.Call("bam_reader_GetReferenceData",object@reader)})
setGeneric("GetReferenceID",function(object,...) standardGeneric("GetReferenceID"))
setMethod(f="GetReferenceID",signature="bamReader",definition=function(object,value) {.Call("bam_reader_GetReferenceID",object@reader,value)})
# This function doesn't work in BamTools
# setGeneric("Jump",function(object,...) standardGeneric("Jump"))
# setMethod(f="Jump",signature="bamReader",definition=function(object,refId,position) { return(.Call("bam_reader_Jump",object@reader,refId,position))})
setGeneric("Rewind",function(object) standardGeneric("Rewind"))
setMethod(f="Rewind",signature="bamReader",definition=function(object) {.Call("bam_reader_Rewind",object@reader)})
setGeneric("CreateIndex",function(object) standardGeneric("CreateIndex"))
setMethod(f="CreateIndex",signature="bamReader",definition=function(object) { return(.Call("bam_reader_CreateIndex",object@reader))})
setGeneric("GetNextAlignment",function(object) standardGeneric("GetNextAlignment"))
setMethod(f="GetNextAlignment",signature="bamReader",definition=function(object) {
							align<-.Call("bam_reader_GetNextAlignment",object@reader)
							align@reader<-object@reader
							return(align)
							})

# #####################################################################################################################
# bamWriter

bamWriter<-function(filename,refData,samHeader) { return(new("bamWriter",filename=filename,refData=refData,samHeader=samHeader))}
setClass("bamWriter", representation(filename="character",writer="externalptr",ref_data="data.frame",sam_header="character"),
			validity=function(object){
				return(TRUE)
			}
)

setMethod(f="initialize", signature="bamWriter",
		definition=function(.Object,filename,refData,samHeader){
			.Object@filename<-filename
			.Object@writer=.Call("get_bam_writer")
			.Object@ref_data<-refData
			.Object@sam_header<-samHeader
			return(.Object)
		}
)

# generic in header
setMethod("filename", "bamWriter", function(object) return(object@filename))
# generic in header
setReplaceMethod(f="filename", signature="bamWriter",
   definition=function(object,value){
       if(!is.character(value))
           stop("class bamWriter, filename setter: value must be character")
       object@filename<-value
       return(object)
   }
)
# generic in header
setMethod(f="bamOpen",signature="bamWriter", definition=function(object) {
					.Call("bam_writer_open",object@writer,object@filename,object@sam_header,object@ref_data)
					return(invisible())
					})
# generic in header
setMethod(f="bamClose",signature="bamWriter",definition=function(object) {.Call("bam_writer_close",object@writer); return(invisible()) })
# generic in header
setMethod("isOpen",signature="bamWriter",definition=function(object) {.Call("bam_writer_isOpen",object@writer)})

# generic in header
setMethod("samHeader", "bamWriter", function(object) return(object@ref_data))
# generic in header
setReplaceMethod(f="samHeader", signature="bamWriter",
   definition=function(object,value){
       if(!is.character(value))
           stop("class bamWriter, samHeader setter: value must be character")
       object@sam_header<-value
       return(object)
   }
)

setGeneric("refData", function(object) standardGeneric("refData"))
setGeneric("refData<-", function(object,value) standardGeneric("refData<-"))
setMethod("refData", "bamWriter", function(object) return(object@ref_data))
setReplaceMethod(f="refData", signature="bamWriter",
   definition=function(object,value){
       if(!is.data.frame(value))
           stop("class bamWriter, refData setter: value must be data.frame")
       object@ref_data<-value
       return(object)
   }
)

setGeneric("bamSave",function(object,...) standardGeneric("bamSave"))
setMethod(f="bamSave",signature="bamWriter",definition=function(object,value) {
							.Call("bam_writer_SaveAlignment",object@writer,value@bamAlignment)
							return(invisible())
})


# ###############################################################################################################################
# BamAlignment

# valid is set directly from c-code in GetNextAlignment
setClass("bamAlignment", representation(bamAlignment="externalptr",valid="logical",reader="externalptr"),
			validity=function(object){
				return(TRUE)
			}
)
setMethod(f="initialize", signature="bamAlignment",
		definition=function(.Object,alignment){
			.Object@bamAlignment<-alignment
			.Object@valid<-FALSE
			.Object@reader<-NULL
			return(.Object)
		}
)


# BamAlignment Member Reader functions
setGeneric("isValid",function(object) standardGeneric("isValid"))
setMethod(f="isValid",signature="bamAlignment",definition=function(object) {return(object@valid)})
setGeneric("refID",function(object) standardGeneric("refID"))
setMethod(f="refID",signature="bamAlignment",definition=function(object) {return(.Call("bam_alignment_getRefID",object@bamAlignment))})
setGeneric("position",function(object) standardGeneric("position"))
setMethod(f="position",signature="bamAlignment",definition=function(object) {return(.Call("bam_alignment_getPosition",object@bamAlignment))})
setGeneric("CigarData",function(object) standardGeneric("CigarData"))
setMethod(f="CigarData",signature="bamAlignment",definition=function(object) { .Call("bam_alignment_getCigar_df",object@bamAlignment)})
setGeneric("MateRefID",function(object) standardGeneric("MateRefID"))
setMethod(f="MateRefID",signature="bamAlignment",definition=function(object) { .Call("bam_alignment_getMateRefID",object@bamAlignment)})
setGeneric("MatePosition",function(object) standardGeneric("MatePosition"))
setMethod(f="MatePosition",signature="bamAlignment",definition=function(object) { .Call("bam_alignment_getMatePosition",object@bamAlignment)})
setGeneric("InsertSize",function(object) standardGeneric("InsertSize"))
setMethod(f="InsertSize",signature="bamAlignment",definition=function(object) { .Call("bam_alignment_getInsertSize",object@bamAlignment)})
setGeneric("MapQuality",function(object) standardGeneric("MapQuality"))
setMethod(f="MapQuality",signature="bamAlignment",definition=function(object) { .Call("bam_alignment_getMapQuality",object@bamAlignment)})
setGeneric("QueryBases",function(object) standardGeneric("QueryBases"))
setMethod(f="QueryBases",signature="bamAlignment",definition=function(object) {return(.Call("bam_alignment_getQueryBases",object@bamAlignment))})
setGeneric("AlignedBases",function(object) standardGeneric("AlignedBases"))
setMethod(f="AlignedBases",signature="bamAlignment",definition=function(object) {return(.Call("bam_alignment_getAlignedBases",object@bamAlignment))})
setGeneric("Qualities",function(object) standardGeneric("Qualities"))
setMethod(f="Qualities",signature="bamAlignment",definition=function(object) {return(.Call("bam_alignment_getQualities",object@bamAlignment))})

# Reads next alignemnt from BAM-file without creating new instances of bamAlignment
setGeneric("getNext",function(.Object) standardGeneric("getNext"))
setMethod(f="getNext",signature="bamAlignment",definition=function(.Object) {
						if(is.null(.Object@reader))
							stop("bamAlignment.getNext: reader slot is uninitialized (-> use GetNextAlignment(reader))\n")
						val<-paste(deparse(substitute(.Object)),"@valid",sep="")
						assign(val,.Call("bam_alignment_getNext",.Object@bamAlignment,.Object@reader),envir=parent.frame())
						return(invisible())
						})


# ##############################################################
# BamAlignment Queries against alignment flag (Readers and Accessors)

setGeneric("name",function(object) standardGeneric("name"))
setMethod(f="name",signature="bamAlignment",definition=function(object) { .Call("bam_alignment_getName",object@bamAlignment) })
setGeneric("name<-",function(object,value) standardGeneric("name<-"))
setReplaceMethod(f="name", signature="bamAlignment",
   definition=function(object,value){
       if(!is.character(value))
           stop("class bamReader, name setter: value must be character")
       .Call("bam_alignment_setName",object@bamAlignment,value)
       return(object)
   }
)


setGeneric("Duplicate", function(object) standardGeneric("Duplicate"))
setGeneric("Duplicate<-", function(object,value) standardGeneric("Duplicate<-"))
setMethod("Duplicate", "bamAlignment", function(object) return(.Call("bam_alignment_IsDuplicate",object@bamAlignment)))
setReplaceMethod(f="Duplicate", signature="bamAlignment",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, Duplicate setter: value must be boolean")
       .Call("bam_alignment_SetIsDuplicate",object@bamAlignment,value)
       return(object)
   }
)

setGeneric("FailedQC", function(object) standardGeneric("FailedQC"))
setGeneric("FailedQC<-", function(object,value) standardGeneric("FailedQC<-"))
setMethod("FailedQC", "bamAlignment", function(object) return(.Call("bam_alignment_IsFailedQC",object@bamAlignment)))
setReplaceMethod(f="FailedQC", signature="bamAlignment",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, FailedQC setter: value must be boolean")
       .Call("bam_alignment_SetIsFailedQC",object@bamAlignment,value)
       return(object)
   }
)

setGeneric("FirstMate", function(object) standardGeneric("FirstMate"))
setGeneric("FirstMate<-", function(object,value) standardGeneric("FirstMate<-"))
setMethod("FirstMate", "bamAlignment", function(object) return(.Call("bam_alignment_IsFirstMate",object@bamAlignment)))
setReplaceMethod(f="FirstMate", signature="bamAlignment",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, FirstMate setter: value must be boolean")
       .Call("bam_alignment_SetIsFirstMate",object@bamAlignment,value)
       return(object)
   }
)

setGeneric("Mapped", function(object) standardGeneric("Mapped"))
setGeneric("Mapped<-", function(object,value) standardGeneric("Mapped<-"))
setMethod("Mapped", "bamAlignment", function(object) return(.Call("bam_alignment_IsMapped",object@bamAlignment)))
setReplaceMethod(f="Mapped", signature="bamAlignment",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, FirstMate setter: value must be boolean")
       .Call("bam_alignment_SetIsUnmapped",object@bamAlignment,(!value))
       return(object)
   }
)

setGeneric("MateMapped", function(object) standardGeneric("MateMapped"))
setGeneric("MateMapped<-", function(object,value) standardGeneric("MateMapped<-"))
setMethod("MateMapped", "bamAlignment", function(object) return(.Call("bam_alignment_IsMateMapped",object@bamAlignment)))
setReplaceMethod(f="MateMapped", signature="bamAlignment",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, FirstMate setter: value must be boolean")
       .Call("bam_alignment_SetIsMateUnmapped",object@bamAlignment,(!value))
       return(object)
   }
)

setGeneric("MateReverseStrand", function(object) standardGeneric("MateReverseStrand"))
setGeneric("MateReverseStrand<-", function(object,value) standardGeneric("MateReverseStrand<-"))
setMethod("MateReverseStrand", "bamAlignment", function(object) return(.Call("bam_alignment_IsMateReverseStrand",object@bamAlignment)))
setReplaceMethod(f="MateReverseStrand", signature="bamAlignment",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, MateReverseStrand setter: value must be boolean")
       .Call("bam_alignment_SetIsMateReverseStrand",object@bamAlignment,value)
       return(object)
   }
)

setGeneric("Paired", function(object) standardGeneric("Paired"))
setGeneric("Paired<-", function(object,value) standardGeneric("Paired<-"))
setMethod("Paired", "bamAlignment", function(object) return(.Call("bam_alignment_IsPaired",object@bamAlignment)))
setReplaceMethod(f="Paired", signature="bamAlignment",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, Paired setter: value must be boolean")
       .Call("bam_alignment_SetIsPaired",object@bamAlignment,value)
       return(object)
   }
)

setGeneric("PrimaryAlignment", function(object) standardGeneric("PrimaryAlignment"))
setGeneric("PrimaryAlignment<-", function(object,value) standardGeneric("PrimaryAlignment<-"))
setMethod("PrimaryAlignment", "bamAlignment", function(object) return(.Call("bam_alignment_IsPrimaryAlignment",object@bamAlignment)))
setReplaceMethod(f="PrimaryAlignment", signature="bamAlignment",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, PrimaryAlignment setter: value must be boolean")
       .Call("bam_alignment_SetIsSecondaryAlignment",object@bamAlignment,(!value))
       return(object)
   }
)

setGeneric("ProperPair", function(object) standardGeneric("ProperPair"))
setGeneric("ProperPair<-", function(object,value) standardGeneric("ProperPair<-"))
setMethod("ProperPair", "bamAlignment", function(object) return(.Call("bam_alignment_IsProperPair",object@bamAlignment)))
setReplaceMethod(f="ProperPair", signature="bamAlignment",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, ProperPair setter: value must be boolean")
       .Call("bam_alignment_SetIsProperPair",object@bamAlignment,value)
       return(object)
   }
)

setGeneric("ReverseStrand", function(object) standardGeneric("ReverseStrand"))
setGeneric("ReverseStrand<-", function(object,value) standardGeneric("ReverseStrand<-"))
setMethod("ReverseStrand", "bamAlignment", function(object) return(.Call("bam_alignment_IsReverseStrand",object@bamAlignment)))
setReplaceMethod(f="ReverseStrand", signature="bamAlignment",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, ReverseStrand setter: value must be boolean")
       .Call("bam_alignment_SetIsReverseStrand",object@bamAlignment,value)
       return(object)
   }
)

setGeneric("SecondMate", function(object) standardGeneric("SecondMate"))
setGeneric("SecondMate<-", function(object,value) standardGeneric("SecondMate<-"))
setMethod("SecondMate", "bamAlignment", function(object) return(.Call("bam_alignment_IsSecondMate",object@bamAlignment)))
setReplaceMethod(f="SecondMate", signature="bamAlignment",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, SecondMate setter: value must be boolean")
       .Call("bam_alignment_SetIsSecondMate",object@bamAlignment,value)
       return(object)
   }
)

