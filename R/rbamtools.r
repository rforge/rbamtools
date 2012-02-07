
# .First.lib(libname, pkgname) is called when a package is loaded; .Last.lib(libpath) is called when a package is detached.

.First.lib<-function(libname,pkgname) { library.dynam("rbamtools", pkgname, libname) }
.Last.lib<-function(libpath) { library.dynam.unload("rbamtools",libpath) }

# #####################################################################################################################
# Generics definitions for bamReader and bamWriter
setGeneric("filename", function(object) standardGeneric("filename"))
setGeneric("isOpen",function(con,rw="") standardGeneric("isOpen"))
setGeneric("bamClose", function(object) standardGeneric("bamClose"))
setGeneric("getNextAlign",function(object) standardGeneric("getNextAlign")) # generic for bamReader and bamRange
setGeneric("bamSave",function(object,...) standardGeneric("bamSave"))       # generif for bamWriter and bamRange

bamReader<-function(filename) { return(new("bamReader",filename)) }

# #####################################################################################################################
# bamReader
setClass("bamReader", representation(filename="character",reader="externalptr",index="externalptr"),
			validity=function(object){
				return(TRUE)
			}
)

setMethod(f="initialize", signature="bamReader",
		definition=function(.Object,filename){
			.Object@filename<-filename
			.Object@reader=.Call("bam_reader_open",filename)
			return(.Object)
		}
)

# generic in header
setMethod("filename", "bamReader", function(object) return(object@filename))
#setGeneric("isOpen",function(con,rw="") standardGeneric("isOpen"))
setMethod("isOpen",signature="bamReader",definition=function(con,rw="") { return(!(.Call("is_nil_externalptr",con@reader,PACKAGE="rbamtools"))) })


setGeneric("index.initialized",function(object) standardGeneric("index.initialized"))
setMethod("index.initialized", signature="bamReader",definition=function(object) {return(!(.Call("is_nil_externalptr",object@index,PACKAGE="rbamtools")))})

# generic in header
setMethod(f="bamClose",signature="bamReader",definition=function(object) {
  .Call("bam_reader_unload_index",object@index,PACKAGE="rbamtools")
  invisible(.Call("bam_reader_close",object@reader,PACKAGE="rbamtools"))
  })
setGeneric("getHeaderText", function(object) standardGeneric("getHeaderText"))
setMethod(f="getHeaderText",signature="bamReader",definition=function(object) { .Call("bam_reader_get_header_text",object@reader,PACKAGE="rbamtools") })
setGeneric("getRefCount",function(object) standardGeneric("getRefCount"))
setMethod(f="getRefCount",signature="bamReader",definition=function(object) {return(.Call("bam_reader_get_ref_count",object@reader,PACKAGE="rbamtools"))})
setGeneric("getRefData",function(object) standardGeneric("getRefData"))
setMethod(f="getRefData",signature="bamReader",definition=function(object) {return(.Call("bam_reader_get_ref_data",object@reader,PACKAGE="rbamtools"))})
setGeneric("createIndex",function(object,idx_filename) standardGeneric("createIndex"))
setMethod(f="createIndex",signature="bamReader",definition=function(object,idx_filename) {invisible(.Call("bam_reader_create_index",object@filename,idx_filename,PACKAGE="rbamtools"))})
# setGeneric("getNextAlign",function(object) standardGeneric("getNextAlign"))
setMethod(f="getNextAlign",signature="bamReader",definition=function(object) { return(new("bamAlign",.Call("bam_reader_get_next_align",object@reader,PACKAGE="rbamtools"))) })

setGeneric("index",function(object) standardGeneric("index"))
setMethod(f="index",signature="bamReader",definition=function(object) {return(!(.Call("is_nil_externalptr",object@index,PACKAGE="rbamtools")))})
setGeneric("index<-",function(object,value) standardGeneric("index<-"))
setReplaceMethod(f="index",signature="bamReader", 
                 definition=function(object,value) {
                   if(!is.character(value))
                     stop("bamReader set index: value must be character!\n")
                   if(!file.exists(value))
                      stop("bamReader set index: index file does not exsist!\n")
                   object@index<-.Call("bam_reader_load_index",value,PACKAGE="rbamtools")
                   return(object)
                   })

setGeneric("bamSort",function(object,prefix,byName=FALSE,maxmem=500000000) standardGeneric("bamSort"))
setMethod(f="bamSort",signature="bamReader",definition=function(object,prefix="sorted",byName=FALSE,maxmem=500000000)
  {
      if(!is.logical(byName))
        stop("[bamSort] byName must be logical!\n")
      if(!is.numeric(maxmem))
        stop("[bamSort] maxmem must be numeric!\n")
      
      cat("[bamSort] Filename: ",object@filename,"\n")
      cat("[bamSort] Prefix:   ",prefix,"\n")
      
      if(byName)
        { cat("[bamSort] Sorting by Name\n") }
      else
        { cat("[bamSort] Sorting by (refID,position)\n") }
      
      cat("[bamSort] MaxMem:   ",maxmem,"\n")
      .Call("bam_reader_sort_file",object@filename,prefix,as.integer(maxmem),byName,PACKAGE="rbamtools")
      cat("[bamSort] Sorting finished.\n")
      invisible(NULL)
  })

# #####################################################################################################################
# bamWriter
setClass("bamWriter", representation(filename="character",writer="externalptr"),
			validity=function(object){
				return(TRUE)
			}
)

setMethod(f="initialize", signature="bamWriter",
		definition=function(.Object,reader,filename){
			.Object@filename<-filename
      .Object@writer<-.Call("bam_writer_open",reader@reader,filename,PACKAGE="rbamtools")
			return(.Object)
		}
)

bamWriter<-function(reader,filename) { return(new("bamWriter",reader,filename)) }

# generic in header
setMethod("filename", "bamWriter", function(object) return(object@filename))
#setGeneric("isOpen",function(con,rw="") standardGeneric("isOpen"))
setMethod("isOpen",signature="bamWriter",definition=function(con,rw="") { return(!(.Call("is_nil_externalptr",con@writer,PACKAGE="rbamtools"))) })


# generic in header
setMethod(f="bamClose",signature="bamWriter",definition=function(object) { invisible(.Call("bam_writer_close",object@writer,PACKAGE="rbamtools"))})
# setGeneric("bamSave",function(object,...) standardGeneric("bamSave"))
setMethod(f="bamSave",signature="bamWriter",definition=function(object,value) 
  {
    if(is(value,"bamAlign"))
    {
      return(invisible(.Call("bam_writer_save_align",object@writer,value@align,PACKAGE="rbamtools")))
    }
    if(is(value,"bamRange"))
    {
      return(invisible(.Call("bam_range_write",object@writer,value@range,PACKAGE="rbamtools")))
    }
    else
      stop("bamSave: Saved object must be of type bamAlign or bamRange!\n")

  })

###################################################################################################
# bamRange

bamRange<-function(reader,coords) { return(new("bamRange",reader,coords))}

setClass("bamRange",representation(range="externalptr"),
         validity=function(object) { return(ifelse(is.null(.Object@range),FALSE,TRUE)) })

setMethod(f="initialize",signature="bamRange",
          definition=function(.Object,reader,coords) {
              if(!is(reader,"bamReader"))
                stop("bamRange initialize: reader must be an instance of bamReader!\n")
              if(!is.numeric(coords)|length(coords)!=3)
                stop("bamRange initialize: coords must be 3-dim integer (ref,start,stop)!\n")
              
              if(is.null(reader@index))
                stop("bamRange initialize: bamReader must have initialized index!\n")
              .Object@range<-.Call("bam_range_fetch",reader@reader,reader@index,as.integer(coords),PACKAGE="rbamtools")
              return(.Object)
              })
setGeneric("as.data.frame",function(x,row.names=NULL,optional=FALSE,...) standardGeneric("as.data.frame"))
setMethod("as.data.frame",signature="bamRange",definition=function(x,row.names=NULL,optional=FALSE,...) {
        df<-.Call("bam_range_get_align_df",x@range,PACKAGE="rbamtools")
        as.data.frame(df,row.names=row.names(df),optional=optional,...)
        })

setGeneric("size",function(object) standardGeneric("size"))
setMethod("size",signature="bamRange",definition=function(object) { .Call("bam_range_get_size",object@range,PACKAGE="rbamtools")} )

#setGeneric("getNextAlign",function(object) standardGeneric("getNextAlign"))
setMethod("getNextAlign",signature="bamRange",definition=function(object) { return(new("bamAlign",.Call("bam_range_get_next_align",object@range,PACKAGE="rbamtools")))})

setGeneric("getPrevAlign",function(object) standardGeneric("getPrevAlign"))
setMethod("getPrevAlign",signature="bamRange",definition=function(object) { return(new("bamAlign",.Call("bam_range_get_prev_align",object@range,PACKAGE="rbamtools")))})

setGeneric("windBack", function(object) standardGeneric("windBack"))
setMethod("windBack",signature="bamRange",definition=function(object) {invisible(.Call("bam_range_wind_back",object@range,PACKAGE="rbamtools"))})

setGeneric("push_back",function(object,value) standardGeneric("push_back"))
setMethod("push_back",signature="bamRange",definition=function(object,value)
  {
    if(!is(value,"bamAlign"))
      stop("push_back.bamRange: pushed object must be of class \"bamAlign\"\n")
    .Call("bam_range_push_back",object@range,value@align,PACKAGE="rbamtools")
  })

setGeneric("pop_back",function(object) standardGeneric("pop_back"))
setMethod("pop_back",signature="bamRange",definition=function(object){.Call("bam_range_pop_back",object@range,PACKAGE="rbamtools") })

setGeneric("push_front",function(object,value) standardGeneric("push_front"))
setMethod("push_front",signature="bamRange",definition=function(object,value)
  {
    if(!is(value,"bamAlign"))
      stop("push_front.bamRange: pushed object must be of class \"bamAlign\"\n")
    .Call("bam_range_push_front",object@range,value@align,PACKAGE="rbamtools")
  })

setGeneric("pop_front",function(object) standardGeneric("pop_front"))
setMethod("pop_front",signature="bamRange",definition=function(object) {.Call("bam_range_pop_front",object@range,PACKAGE="rbamtools")})

setGeneric("writeCurrentAlign",function(object,value) standardGeneric("writeCurrentAlign"))
setMethod("writeCurrentAlign",signature="bamRange",definition=function(object,value)
  {
      if(!is(value,"bamAlign"))
        stop("writeCurrentAlign.bamRange: written object must be of class \"bamAlign\"\n")
      .Call("bam_range_write_current_align",object@range,value@align,PACKAGE="rbamtools")
  })

setGeneric("insertPastCurrent",function(object,value) standardGeneric("insertPastCurrent"))
setMethod("insertPastCurrent",signature="bamRange",definition=function(object,value)
  {
      if(!is(value,"bamAlign"))
        stop("insertPastCurrent.bamRange: written object must be of class \"bamAlign\"\n")
      .Call("bam_range_insert_past_curr_align",object@range,value@align,PACKAGE="rbamtools")
  })

setGeneric("insertPreCurrent",function(object,value) standardGeneric("insertPreCurrent"))
setMethod("insertPreCurrent",signature="bamRange",definition=function(object,value)
  {
      if(!is(value,"bamAlign"))
        stop("insertPreCurrent.bamRange: written object must be of class \"bamAlign\"\n")
      .Call("bam_range_insert_pre_curr_align",object@range,value@align,PACKAGE="rbamtools")
  })



# ###############################################################################################################################
# bamAlign
setClass("bamAlign", representation(align="externalptr"),
			validity=function(object){
				return(ifelse(is.null(.Object@align,FALSE,TRUE)))
			}
)
setMethod(f="initialize", signature="bamAlign",
		definition=function(.Object,alignment){
			.Object@align<-alignment
			return(.Object)
		}
)

# bamAlign Member Reader functions
setGeneric("name",function(object) standardGeneric("name"))
setMethod(f="name",signature="bamAlign",definition=function(object) { .Call("bam_alignment_get_name",object@align,PACKAGE="rbamtools") })
setGeneric("refID",function(object) standardGeneric("refID"))
setMethod(f="refID",signature="bamAlign",definition=function(object) {return(.Call("bam_alignment_get_refid",object@align,PACKAGE="rbamtools"))})
setGeneric("position",function(object) standardGeneric("position"))
setMethod(f="position",signature="bamAlign",definition=function(object) {return(.Call("bam_alignment_get_position",object@align,PACKAGE="rbamtools"))})
setGeneric("cigarData",function(object) standardGeneric("cigarData"))
setMethod(f="cigarData",signature="bamAlign",definition=function(object) { .Call("bam_alignment_get_cigar_df",object@align,PACKAGE="rbamtools")})
setGeneric("mateRefID",function(object) standardGeneric("mateRefID"))
setMethod(f="mateRefID",signature="bamAlign",definition=function(object) { .Call("bam_alignment_get_mate_refid",object@align,PACKAGE="rbamtools")})
setGeneric("matePosition",function(object) standardGeneric("matePosition"))
setMethod(f="matePosition",signature="bamAlign",definition=function(object) { .Call("bam_alignment_get_mate_position",object@align,PACKAGE="rbamtools")})
setGeneric("insertSize",function(object) standardGeneric("insertSize"))
setMethod(f="insertSize",signature="bamAlign",definition=function(object) { .Call("bam_alignment_get_insert_size",object@align,PACKAGE="rbamtools")})
setGeneric("mapQuality",function(object) standardGeneric("mapQuality"))
setMethod(f="mapQuality",signature="bamAlign",definition=function(object) { .Call("bam_alignment_get_map_quality",object@align,PACKAGE="rbamtools")})
setGeneric("readBases",function(object) standardGeneric("readBases"))
setMethod(f="readBases",signature="bamAlign",definition=function(object) {return(.Call("bam_alignment_get_read_bases",object@align,PACKAGE="rbamtools"))})
setGeneric("qualities",function(object) standardGeneric("qualities"))
setMethod(f="qualities",signature="bamAlign",definition=function(object) {return(.Call("bam_alignment_get_qualities",object@align,PACKAGE="rbamtools"))})

# ##############################################################
# Queries against alignment flag (Readers and Accessors)

setGeneric("pcrORopt_duplicate", function(object) standardGeneric("pcrORopt_duplicate"))
setGeneric("pcrORopt_duplicate<-", function(object,value) standardGeneric("pcrORopt_duplicate<-"))
setMethod("pcrORopt_duplicate", "bamAlign", function(object) return(.Call("bam_alignment_is_pcr_or_optical_dup",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="pcrORopt_duplicate", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, Duplicate setter: value must be boolean")
         .Call("bam_alignment_set_is_pcr_or_optical_dup",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("failedQC", function(object) standardGeneric("failedQC"))
setGeneric("failedQC<-", function(object,value) standardGeneric("failedQC<-"))
setMethod("failedQC", "bamAlign", function(object) return(.Call("bam_alignment_fail_qc",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="failedQC", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, failedQC setter: value must be boolean")
       .Call("bam_alignment_set_fail_qc",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("firstInPair", function(object) standardGeneric("firstInPair"))
setGeneric("firstInPair<-", function(object,value) standardGeneric("firstInPair<-"))
setMethod("firstInPair", "bamAlign", function(object) return(.Call("bam_alignment_is_first_in_pair",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="firstInPair", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, FirstInPair setter: value must be boolean")
       .Call("bam_alignment_set_is_first_in_pair",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("secondInPair", function(object) standardGeneric("secondInPair"))
setGeneric("secondInPair<-", function(object,value) standardGeneric("secondInPair<-"))
setMethod("secondInPair", "bamAlign", function(object) return(.Call("bam_alignment_is_second_in_pair",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="secondInPair", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, secondInPair setter: value must be boolean")
       .Call("bam_alignment_set_is_second_in_pair",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("unmapped", function(object) standardGeneric("unmapped"))
setGeneric("unmapped<-", function(object,value) standardGeneric("unmapped<-"))
setMethod("unmapped", "bamAlign", function(object) return(.Call("bam_alignment_is_unmapped",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="unmapped", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, unmapped setter: value must be boolean")
       .Call("bam_alignment_set_is_unmapped",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("mateUnmapped", function(object) standardGeneric("mateUnmapped"))
setGeneric("mateUnmapped<-", function(object,value) standardGeneric("mateUnmapped<-"))
setMethod("mateUnmapped", "bamAlign", function(object) return(.Call("bam_alignment_mate_is_unmapped",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="mateUnmapped", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, mateUnmapped setter: value must be boolean")
       .Call("bam_alignment_set_mate_is_unmapped",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("reverseStrand", function(object) standardGeneric("reverseStrand"))
setGeneric("reverseStrand<-", function(object,value) standardGeneric("reverseStrand<-"))
setMethod("reverseStrand", "bamAlign", function(object) return(.Call("bam_alignment_strand_reverse",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="reverseStrand", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, reverseStrand setter: value must be boolean")
       .Call("bam_alignment_set_strand_reverse",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("mateReverseStrand", function(object) standardGeneric("mateReverseStrand"))
setGeneric("mateReverseStrand<-", function(object,value) standardGeneric("mateReverseStrand<-"))
setMethod("mateReverseStrand", "bamAlign", function(object) return(.Call("bam_alignment_mate_strand_reverse",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="mateReverseStrand", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, mateReverseStrand setter: value must be boolean")
       .Call("bam_alignment_set_mate_strand_reverse",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("paired", function(object) standardGeneric("paired"))
setGeneric("paired<-", function(object,value) standardGeneric("paired<-"))
setMethod("paired", "bamAlign", function(object) return(.Call("bam_alignment_is_paired",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="paired", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, paired setter: value must be boolean")
       .Call("bam_alignment_set_is_paired",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)


setGeneric("properPair", function(object) standardGeneric("properPair"))
setGeneric("properPair<-", function(object,value) standardGeneric("properPair<-"))
setMethod("properPair", "bamAlign", function(object) return(.Call("bam_alignment_mapped_in_proper_pair",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="properPair", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, properPair setter: value must be boolean")
       .Call("bam_alignment_set_mapped_in_proper_pair",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("secondaryAlign", function(object) standardGeneric("secondaryAlign"))
setGeneric("secondaryAlign<-", function(object,value) standardGeneric("secondaryAlign<-"))
setMethod("secondaryAlign", "bamAlign", function(object) return(.Call("bam_alignment_is_secondary_align",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="secondaryAlign", signature="bamAlign",
   definition=function(object,value){
       if(!is.logical(value))
           stop("class bamReader, SecondaryAlign setter: value must be boolean")
       .Call("bam_alignment_set_is_secondary_align",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)

setGeneric("flag", function(object) standardGeneric("flag"))
setGeneric("flag<-", function(object,value) standardGeneric("flag<-"))
setMethod("flag", "bamAlign", function(object) return(.Call("bam_alignment_get_flag",object@align,PACKAGE="rbamtools")))
setReplaceMethod(f="flag", signature="bamAlign",
   definition=function(object,value){
       if(!is.integer(value))
           stop("class bamReader, flag setter: value must be boolean")
       .Call("bam_alignment_set_flag",object@align,value,PACKAGE="rbamtools")
       return(object)
   }
)
