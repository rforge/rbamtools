/*
 * rbamtools.cpp
 *
 *  Created on: 17.06.2011
 *      Author: wolfgang
 */

#include "rbamtools.h"

extern "C" {

	///////////////////////////////////////////////////////////////////////////////////////////////
	// BamReader
	///////////////////////////////////////////////////////////////////////////////////////////////
	static void finalize_bam_reader(SEXP ptr)
	{
		if(TYPEOF(ptr)!=EXTPTRSXP)
		{
			error("finalize_bam_reader: no external pointer!");
			return;
		}
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(ptr));
		if(reader)
		{
			reader->Close();
			delete reader;
			reader=0;
			R_SetExternalPtrAddr(ptr,NULL);
			Rprintf("finalize_bam_reader: BamReader finalized\n");
		}
		else
		{
			Rprintf("finalize_bam_reader: Nothing to finalize\n");
		}
	}

	SEXP get_bam_reader()
	{
		BamReader *reader=new BamReader();
		SEXP ptr;
		PROTECT(ptr=R_MakeExternalPtr(static_cast<void*>(reader),R_NilValue,R_NilValue));
		R_RegisterCFinalizer(ptr,finalize_bam_reader);
		UNPROTECT(1);
		return ptr;
	}

	SEXP bam_reader_open(SEXP pReader,SEXP filename,SEXP index_file_name)
	{
		if(TYPEOF(pReader)!=EXTPTRSXP)
		{
			error("bam_connect: no external pointer!");
			return R_NilValue;
		}
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));

		if(TYPEOF(filename)!=STRSXP)
		{
			error("bam_connect: filename must be a string");
			return R_NilValue;
		}

		if(TYPEOF(index_file_name)!=STRSXP)
		{
			error("bam_connect: index_file_name must be a string");
			return R_NilValue;
		}

		const char *c=CHAR(STRING_ELT(filename,0));
		const char *i=CHAR(STRING_ELT(index_file_name,0));
		reader->Open(c,i);
		Rprintf("BamReader opened file %s\n",c);
		return R_NilValue;
	}

	SEXP bam_reader_isOpen(SEXP pReader)
	{
		if(TYPEOF(pReader)!=EXTPTRSXP)
		{
			error("bam_reader_isOpen: no external pointer!\n");
		}
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));

		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=reader->isOpen();
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_reader_close(SEXP pReader)
	{
		if(TYPEOF(pReader)!=EXTPTRSXP)
		{
			error("bam_disconnect: no external pointer!");
			return R_NilValue;
		}
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));
		reader->Close();
		Rprintf("BamReader closed\n");
		return R_NilValue;
	}

	SEXP bam_reader_GetHeaderText(SEXP pReader)
	{
		if(TYPEOF(pReader)!=EXTPTRSXP)
		{
			error("bam_reader_GetHeaderText: no external pointer!");
			return R_NilValue;
		}
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));
		if(!(reader->isOpen()))
		{
			error("bam_reader_GetHeaderText: BAM file is not open!\n");
		}

		SEXP ans;
		PROTECT(ans=Rf_allocVector(STRSXP,1));
		SET_STRING_ELT(ans,0,mkChar(reader->GetHeaderText().c_str()));
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_reader_GetReferenceCount(SEXP pReader)
	{
		if(TYPEOF(pReader)!=EXTPTRSXP)
		{
			error("bam_reader_GetReferenceCount: no external pointer!");
			return R_NilValue;
		}
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));
		if(!(reader->isOpen()))
		{
			error("bam_reader_GetReferenceCount: BAM file is not open!\n");
		}
		SEXP ans;
		PROTECT(ans=Rf_allocVector(INTSXP,1));
		INTEGER(ans)[0]=reader->GetReferenceCount();
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_reader_GetReferenceData(SEXP pReader)
	{
		if(TYPEOF(pReader)!=EXTPTRSXP)
		{
			error("bam_reader_GetReferenceData: no external pointer!");
			return R_NilValue;
		}
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));
		if(!(reader->isOpen()))
		{
			error("bam_reader_GetReferenceData: BAM file is not open!\n");
		}
		vector<RefData> v=reader->GetReferenceData();

		// create data.frame
		int nProtected=0;
		int nCols=4;
		SEXP dflist;
		PROTECT(dflist=allocVector(VECSXP,nCols));
		++nProtected;
		int nRows=v.size();

		// Column 0: RefID
		SEXP RefID_vector;
		PROTECT(RefID_vector=allocVector(INTSXP,nRows));
		++nProtected;

		// Column 1: RefName
		SEXP RefName_vector;
		PROTECT(RefName_vector=allocVector(STRSXP,nRows));
		++nProtected;

		// Column 2: RefLength
		SEXP RefLength_vector;
		PROTECT(RefLength_vector=allocVector(INTSXP,nRows));
		++nProtected;

		// Column 3: RefHasAlignments
		SEXP RefHas_vector;
		PROTECT(RefHas_vector=allocVector(LGLSXP,nRows));
		++nProtected;

		for(int i=0;i<nRows;++i)
		{
			INTEGER(RefID_vector)[i]=i;
			SET_STRING_ELT(RefName_vector,i,mkChar(v[i].RefName.c_str()));
			INTEGER(RefLength_vector)[i]=v[i].RefLength;
			LOGICAL(RefHas_vector)[i]=v[i].RefHasAlignments;
		}
		SET_VECTOR_ELT(dflist,0,RefID_vector);
		SET_VECTOR_ELT(dflist,1,RefName_vector);
		SET_VECTOR_ELT(dflist,2,RefLength_vector);
		SET_VECTOR_ELT(dflist,3,RefHas_vector);

		// Column Names
		SEXP col_names;
		PROTECT(col_names=allocVector(STRSXP,nCols));
		++nProtected;
		SET_STRING_ELT(col_names,0,mkChar("RefID"));
		SET_STRING_ELT(col_names,1,mkChar("RefName"));
		SET_STRING_ELT(col_names,2,mkChar("RefLength"));
		SET_STRING_ELT(col_names,3,mkChar("RefHasAlignments"));
		setAttrib(dflist,R_NamesSymbol,col_names);

		SEXP row_names;
	    PROTECT(row_names=allocVector(STRSXP,nRows));
	    ++nProtected;
	    char c[20];
	    for(int i=0;i<nRows;++i)
	    {
	    	sprintf(c,"%i",i);
	    	SET_STRING_ELT(row_names,i,mkChar(c));
	    }
        setAttrib(dflist,R_RowNamesSymbol,row_names);

		setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
		UNPROTECT(nProtected);
		return dflist;
	}

	SEXP bam_reader_GetReferenceID(SEXP pReader, SEXP refName)
	{
		if(TYPEOF(pReader)!=EXTPTRSXP)
		{
			error("bam_reader_GetReferenceID: no external pointer!");
			return R_NilValue;
		}
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));
		if(!(reader->isOpen()))
		{
			error("bam_reader_GetReferenceID: BAM file is not open!\n");
		}

		// Read search string from argument
		string sRefName(CHAR(STRING_ELT(refName,0)));
		SEXP ans;
		PROTECT(ans=Rf_allocVector(INTSXP,1));
		INTEGER(ans)[0]=reader->GetReferenceID(sRefName);
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_reader_Jump(SEXP pReader,SEXP refID,SEXP position)
	{
		if(TYPEOF(pReader)!=EXTPTRSXP)
		{
			error("bam_reader_Jump: no external pointer!");
			return R_NilValue;
		}
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));
		if(!(reader->isOpen()))
		{
			error("bam_reader_Jump: BAM file is not open!\n");
		}

		int *pRefID=INTEGER_POINTER(AS_INTEGER(refID));
		int *pPosition=INTEGER_POINTER(AS_INTEGER(position));

        	SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));

		// bamtools
		LOGICAL(ans)[0]=reader->Jump(*pRefID,*pPosition);
        	UNPROTECT(1);
       		return ans;
	}

	SEXP bam_reader_Rewind(SEXP pReader)
	{
		if(TYPEOF(pReader)!=EXTPTRSXP)
		{
			error("bam_reader_Rewind: no external pointer!");
			return R_NilValue;
		}
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));
		if(!(reader->isOpen()))
		{
			error("bam_reader_Rewind: BAM file is not open\n");
		}

		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		// bamtools
		LOGICAL(ans)[0]=reader->Rewind();
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_reader_CreateIndex(SEXP pReader)
	{
		if(TYPEOF(pReader)!=EXTPTRSXP)
		{
			error("bam_reader_CreateIndex: no external pointer!\n");
			return R_NilValue;
		}
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));
		if(!(reader->isOpen()))
		{
			error("bam_reader_CreateIndex: BAM file is not open!\n");
		}


		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		// bamtools
		LOGICAL(ans)[0]=reader->CreateIndex();
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_reader_GetNextAlignment(SEXP pReader)
	{
		if(TYPEOF(pReader)!=EXTPTRSXP)
			error("bam_reader_GetNextAlignment: no external pointer!");
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));
		if(!(reader->isOpen()))
		{
			error("bam_reader_GetNextAlignment: BAM file is not open!\n");
		}

		BamAlignment *pAlign=new BamAlignment;
		int nProtected=0;

		SEXP ans,ptr,valid;
		PROTECT(ans=NEW_OBJECT(MAKE_CLASS("bamAlignment")));
		++nProtected;
		PROTECT(valid=Rf_allocVector(LGLSXP,1));
		++nProtected;

		if(reader->GetNextAlignment(*pAlign))
		{
			LOGICAL(valid)[0]=true;
		}
		else
		{
			LOGICAL(valid)[0]=false;
		}
		PROTECT(ptr=R_MakeExternalPtr(static_cast<void*>(pAlign),R_NilValue,R_NilValue));
		++nProtected;
		R_RegisterCFinalizer(ptr,finalize_BamAlignment);

		SEXP ptrSlotName,validSlotName;
		PROTECT(ptrSlotName=Rf_mkString("bamAlignment"));
		++nProtected;
		PROTECT(validSlotName=Rf_mkString("valid"));
		++nProtected;

		ans=SET_SLOT(ans,ptrSlotName,ptr);
		ans=SET_SLOT(ans,validSlotName,valid);

		UNPROTECT(nProtected);
		return ans;
	}


	///////////////////////////////////////////////////////////////////////////////////////////////
	// BamWriter
	///////////////////////////////////////////////////////////////////////////////////////////////

	static void finalize_bam_writer(SEXP ptr)
	{
		if(TYPEOF(ptr)!=EXTPTRSXP)
		{
			error("finalize_bam_writer: no external pointer!");
			return;
		}
		BamWriter *writer=static_cast<BamWriter*>(R_ExternalPtrAddr(ptr));
		if(writer)
		{
			writer->Close();
			delete writer;
			writer=0;
			R_SetExternalPtrAddr(ptr,NULL);
			Rprintf("finalize_bam_writer: finalized!\n");
		}
		else
		{
			Rprintf("finalize_bam_writer: nothing to finalize!\n");
		}
	}

	SEXP get_bam_writer()
	{
		BamWriter *writer=new BamWriter();
		SEXP ptr;
		PROTECT(ptr=R_MakeExternalPtr(static_cast<void*>(writer),R_NilValue,R_NilValue));
		R_RegisterCFinalizer(ptr,finalize_bam_writer);
		UNPROTECT(1);
		return ptr;
	}

	SEXP bam_writer_open_reader(SEXP pWriter, SEXP pReader,SEXP pFilename)
	{
		if(TYPEOF(pWriter)!=EXTPTRSXP)
		{
			error("bam_writer_open_reader: pWriter no external pointer!\n");
			return R_NilValue;
		}
		if(TYPEOF(pReader)!=EXTPTRSXP)
		{
			error("bam_writer_open_reader: pReader no external pointer!\n");
			return R_NilValue;
		}
		if(TYPEOF(pFilename)!=STRSXP)
		{
			error("bam_writer_open_reader: pFilename no string!\n");
			return R_NilValue;
		}
		BamWriter *writer=static_cast<BamWriter*>(R_ExternalPtrAddr(pWriter));
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));
		if(!(reader->isOpen()))
		{
			error("bam_writer_open_reader: Reader BAM file is not open!\n");
			return R_NilValue;
		}
		writer->Open(string(CHAR(STRING_ELT(pFilename,0))),reader->GetHeaderText(),reader->GetReferenceData());
		return R_NilValue;
	}

	vector<RefData> bam_writer_dfToVector(SEXP pRefData)
	{
		if(TYPEOF(pRefData)!=VECSXP)
		{
			error("bam_writer_dfToVector: no VECSXP\n");
		}
		if(strcmp("data.frame",CHAR(STRING_ELT(getAttrib(pRefData,R_ClassSymbol),0)))!=0)
		{
			error("bam_writer_dfToVector: vector is no data.frame!\n");
		}

		int nProtected=0;
		// Expects four columns: RefID, RefName, RefLength, RefHasAlignments
		// RefID is expected to be 0,1,2,... (and not included into header)
		SEXP RefNames,RefLength,RefHas;
		PROTECT(RefNames=VECTOR_ELT(pRefData,1));
		++nProtected;
		PROTECT(RefLength=VECTOR_ELT(pRefData,2));
		++nProtected;
		PROTECT(RefHas=VECTOR_ELT(pRefData,3));
		++nProtected;

		int nRows=LENGTH(RefNames);
		vector<RefData> v;
		v.reserve(nRows);
		for(int i=0;i<nRows;++i)
		{
			RefData r;
			r.RefName=string(CHAR(STRING_ELT(RefNames,i)));
			r.RefLength=INTEGER(AS_INTEGER(RefLength))[i];
			r.RefHasAlignments=LOGICAL(RefHas)[i];
			v.push_back(r);
		}
		UNPROTECT(nProtected);
		return v;
	}

	SEXP bam_writer_open(SEXP pWriter,SEXP pFilename,SEXP pSamHeader,SEXP pRefSeqs)
	{
		if(TYPEOF(pWriter)!=EXTPTRSXP)
		{
			error("bam_writer_connect: pWriter no external pointer!\n");
			return R_NilValue;
		}
		if(TYPEOF(pSamHeader)!=STRSXP)
		{
			error("bam_writer_connect: pSamHeader not String!\n");
			return R_NilValue;
		}
		if(TYPEOF(pRefSeqs)!=VECSXP)
		{
			error("bam_writer_connect: pRefSeqs not vector!\n");
			return R_NilValue;
		}
		if(strcmp("data.frame",CHAR(STRING_ELT(getAttrib(pRefSeqs,R_ClassSymbol),0)))!=0)
		{
			error("bam_writer_dfToVector: vector is no data.frame!\n");
		}
		BamWriter *writer=static_cast<BamWriter*>(R_ExternalPtrAddr(pWriter));
		writer->Open(string(CHAR(STRING_ELT(pFilename,0))),string(CHAR(STRING_ELT(pSamHeader,0))),bam_writer_dfToVector(pRefSeqs));
		Rprintf("BamWriter opened.\n");
		return R_NilValue;
	}

	SEXP bam_writer_isOpen(SEXP pWriter)
	{
		if(TYPEOF(pWriter)!=EXTPTRSXP)
		{
			error("bam_writer_isOpen: no external pointer!\n");
		}
		BamWriter *writer=static_cast<BamWriter*>(R_ExternalPtrAddr(pWriter));

		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=writer->isOpen();
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_writer_SaveAlignment(SEXP pWriter, SEXP pAlignment)
	{
		if(TYPEOF(pWriter)!=EXTPTRSXP)
		{
			error("bam_writer_SaveAlignment writer: no external pointer!\n");
			return R_NilValue;
		}
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bam_writer_SaveAlignment Alignment: no external pointer!\n");
			return R_NilValue;
		}
		BamWriter *writer=static_cast<BamWriter*>(R_ExternalPtrAddr(pWriter));
		if(!(writer->isOpen()))
		{
			error("bam_writer_SaveAlignment: BAM file is not open!\n");
		}
		BamAlignment *align=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		writer->SaveAlignment(*align);
		return R_NilValue;
	}

	SEXP bam_writer_close(SEXP pWriter)
	{
		if(TYPEOF(pWriter)!=EXTPTRSXP)
		{
			error("bam_writer_close: no exteranl pointer!\n");
			return R_NilValue;
		}
		BamWriter *writer=static_cast<BamWriter*>(R_ExternalPtrAddr(pWriter));
		writer->Close();
		Rprintf("BamWriter closed\n");
		return R_NilValue;
	}


	///////////////////////////////////////////////////////////////////////////////////////////////
	// BamAlignment
	///////////////////////////////////////////////////////////////////////////////////////////////

	static void finalize_BamAlignment(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
			error("BamAlignment_getName: no external pointer!");
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		if(pAlign)
		{
			delete pAlign;
			pAlign=0;
			R_SetExternalPtrAddr(pAlignment,NULL);
			//Rprintf("finalize_BamAlignment: finalized\n");
		}
		else
		{
			Rprintf("finalize_BamAlignment: Nothing to finalize\n");
		}
	}

	SEXP bam_alignment_getNext(SEXP pAlign, SEXP pReader)
	{
		if(TYPEOF(pAlign)!=EXTPTRSXP)
			error("bam_alignment_getNext: no external pointer!");
		if(TYPEOF(pReader)!=EXTPTRSXP)
			error("bam_alignment_getNext: no external pointer!");

		BamAlignment *align=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlign));
		BamReader *reader=static_cast<BamReader*>(R_ExternalPtrAddr(pReader));
		if(!(reader->isOpen()))
		{
			error("bam_reader_GetNextAlignment: BAM file is not open!\n");
		}

		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=(reader->GetNextAlignment(*align));
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_getName(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
			error("BamAlignment_getName: no external pointer!");
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=Rf_allocVector(STRSXP,1));
		SET_STRING_ELT(ans,0,mkChar(pAlign->Name.c_str()));
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_setName(SEXP pAlignment,SEXP newname)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
			error("BamAlignment_setName: no external pointer!");
		if(TYPEOF(newname)!=STRSXP)
			error("bamAlignment_setName: new name must be string!");

		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		(pAlign->Name)=string(CHAR(STRING_ELT(newname,0)));
		return R_NilValue;
	}

	SEXP bam_alignment_getRefID(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
			error("BamAlignment_getRefID: no external pointer!");
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=Rf_allocVector(INTSXP,1));
		INTEGER(ans)[0]=pAlign->RefID;
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_getPosition(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
			error("BamAlignment_getPosition: no external pointer!");
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=Rf_allocVector(INTSXP,1));
		INTEGER(ans)[0]=pAlign->Position;
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_getCigar_df(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
			error("BamAlignment_getCigar_df: no external pointer");
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		vector<CigarOp> v=(pAlign->CigarData);

		// create data.frame
		int nProtected=0;
		int nCols=2;
		SEXP dflist;
		PROTECT(dflist=allocVector(VECSXP,nCols));
		++nProtected;
		int nRows=v.size();
		int i;

		// Column 0: Length
		SEXP Length_vector;
		PROTECT(Length_vector=allocVector(INTSXP,nRows));
		++nProtected;

		// Column 1: Type
		SEXP Type_vector;
		PROTECT(Type_vector=allocVector(STRSXP,nRows));
		++nProtected;

		for(i=0;i<nRows;++i)
		{
			INTEGER(Length_vector)[i]=v[i].Length;
			SET_STRING_ELT(Type_vector,i,mkCharLen(&(v[i].Type),1));
		}

		SET_VECTOR_ELT(dflist,0,Length_vector);
		SET_VECTOR_ELT(dflist,1,Type_vector);

		// Column Names
		SEXP col_names;
		PROTECT(col_names=allocVector(STRSXP,nCols));
		++nProtected;
		SET_STRING_ELT(col_names,0,mkChar("Length"));
		SET_STRING_ELT(col_names,1,mkChar("Type"));
		setAttrib(dflist,R_NamesSymbol,col_names);

		SEXP row_names;
	    PROTECT(row_names=allocVector(STRSXP,nRows));
	    ++nProtected;

	    char c[20];
	    for(i=0;i<nRows;++i)
	    {
	    	sprintf(c,"%i",i);
	    	SET_STRING_ELT(row_names,i,mkChar(c));
	    }
	    setAttrib(dflist,R_RowNamesSymbol,row_names);

		setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
		UNPROTECT(nProtected);
		return dflist;
	}

	SEXP bam_alignment_getMateRefID(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
			error("BamAlignment_getMateRefID: no external pointer!");
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=Rf_allocVector(INTSXP,1));
		INTEGER(ans)[0]=pAlign->MateRefID;
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_getMatePosition(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
			error("BamAlignment_getMatePosition: no external pointer!");
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=Rf_allocVector(INTSXP,1));
		INTEGER(ans)[0]=pAlign->MatePosition;
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_getInsertSize(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
			error("BamAlignment_getInsertSize: no external pointer!");
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=Rf_allocVector(INTSXP,1));
		INTEGER(ans)[0]=pAlign->InsertSize;
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_getMapQuality(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
			error("BamAlignment_getMapQuality: no external pointer!");
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=Rf_allocVector(INTSXP,1));
		INTEGER(ans)[0]=pAlign->MapQuality;
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_getQueryBases(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bam_Alignment_getQueryBases: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=allocVector(STRSXP,1));
		SET_STRING_ELT(ans,0,mkChar((pAlign->QueryBases).c_str()));
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_getAlignedBases(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bam_Alignment_getAlignedBases: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=allocVector(STRSXP,1));
		SET_STRING_ELT(ans,0,mkChar((pAlign->AlignedBases).c_str()));
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_getQualities(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bam_Alignment_getQualities: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=allocVector(STRSXP,1));
		SET_STRING_ELT(ans,0,mkChar((pAlign->Qualities).c_str()));
		UNPROTECT(1);
		return ans;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	// BamAux.h Queries against alignment flag

	SEXP bam_alignment_IsDuplicate(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_IsDuplicate: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=pAlign->IsDuplicate();
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_IsFailedQC(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_IsFailedQC: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));

		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=pAlign->IsFailedQC();
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_IsFirstMate(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_IsFirstMate: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=pAlign->IsFirstMate();
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_IsMapped(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_IsMapped: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=pAlign->IsMapped();
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_IsMateMapped(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_IsMateMapped: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=pAlign->IsMateMapped();
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_IsMateReverseStrand(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_IsMateReverseStrand: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=pAlign->IsMateReverseStrand();
		UNPROTECT(1);
		return ans;
	}

	SEXP bam_alignment_IsPaired(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_IsPaired: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=pAlign->IsPaired();
		UNPROTECT(1);
		return ans;
	}
	SEXP bam_alignment_IsPrimaryAlignment(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_IsPrimaryAlignment: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=pAlign->IsPrimaryAlignment();
		UNPROTECT(1);
		return ans;
	}
	SEXP bam_alignment_IsProperPair(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_IsProperPair: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=pAlign->IsProperPair();
		UNPROTECT(1);
		return ans;
	}
	SEXP bam_alignment_IsReverseStrand(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_IsReverseStrand: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=pAlign->IsReverseStrand();
		UNPROTECT(1);
		return ans;
	}
	SEXP bam_alignment_IsSecondMate(SEXP pAlignment)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_IsSecondMate: no external pointer!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		SEXP ans;
		PROTECT(ans=allocVector(LGLSXP,1));
		LOGICAL(ans)[0]=pAlign->IsSecondMate();
		UNPROTECT(1);
		return ans;
	}

	// Writing accessors
	SEXP bam_alignment_SetIsDuplicate(SEXP pAlignment, SEXP bool_ok)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_SetIsDuplicate: no external pointer!");
			return R_NilValue;
		}
		if(TYPEOF(bool_ok)!=LGLSXP)
		{
			error("bamAlignment_SetIsDupicate: no bool value!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		pAlign->SetIsDuplicate(*(LOGICAL(AS_LOGICAL(bool_ok))));
		return R_NilValue;
	}


	SEXP bam_alignment_SetIsFailedQC(SEXP pAlignment, SEXP bool_ok)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_SetIsFailedQC: no external pointer!");
			return R_NilValue;
		}
		if(TYPEOF(bool_ok)!=LGLSXP)
		{
			error("bamAlignment_SetIsFailedQC: no bool value!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		pAlign->SetIsFailedQC(*(LOGICAL(AS_LOGICAL(bool_ok))));
		return R_NilValue;
	}

	SEXP bam_alignment_SetIsFirstMate(SEXP pAlignment, SEXP bool_ok)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_SetIsFirstMate: no external pointer!");
			return R_NilValue;
		}
		if(TYPEOF(bool_ok)!=LGLSXP)
		{
			error("bamAlignment_SetIsFirstMate: no bool value!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		pAlign->SetIsFirstMate(*(LOGICAL(AS_LOGICAL(bool_ok))));
		return R_NilValue;
	}

	SEXP bam_alignment_SetIsMateUnmapped(SEXP pAlignment, SEXP bool_ok)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_SetIsMateUnmapped: no external pointer!");
			return R_NilValue;
		}
		if(TYPEOF(bool_ok)!=LGLSXP)
		{
			error("bamAlignment_SetIsMateUnmapped: no bool value!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		pAlign->SetIsMateUnmapped(*(LOGICAL(AS_LOGICAL(bool_ok))));
		return R_NilValue;
	}

	SEXP bam_alignment_SetIsMateReverseStrand(SEXP pAlignment, SEXP bool_ok)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_SetIsMateReverseStrand: no external pointer!");
			return R_NilValue;
		}
		if(TYPEOF(bool_ok)!=LGLSXP)
		{
			error("bamAlignment_SetIsMateReverseStrand: no bool value!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		pAlign->SetIsMateReverseStrand(*(LOGICAL(AS_LOGICAL(bool_ok))));
		return R_NilValue;
	}

	SEXP bam_alignment_SetIsPaired(SEXP pAlignment, SEXP bool_ok)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_SetIsPaired: no external pointer!");
			return R_NilValue;
		}
		if(TYPEOF(bool_ok)!=LGLSXP)
		{
			error("bamAlignment_SetIsPaired: no bool value!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		pAlign->SetIsPaired(*(LOGICAL(AS_LOGICAL(bool_ok))));
		return R_NilValue;
	}

	SEXP bam_alignment_SetIsProperPair(SEXP pAlignment, SEXP bool_ok)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_SetIsProperPair: no external pointer!");
			return R_NilValue;
		}
		if(TYPEOF(bool_ok)!=LGLSXP)
		{
			error("bamAlignment_SetIsProperPair: no bool value!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		pAlign->SetIsProperPair(*(LOGICAL(AS_LOGICAL(bool_ok))));
		return R_NilValue;
	}
	//
	SEXP bam_alignment_SetIsReverseStrand(SEXP pAlignment, SEXP bool_ok)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_SetIsReverseStrand: no external pointer!");
			return R_NilValue;
		}
		if(TYPEOF(bool_ok)!=LGLSXP)
		{
			error("bamAlignment_SetIsReverseStrand: no bool value!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		pAlign->SetIsReverseStrand(*(LOGICAL(AS_LOGICAL(bool_ok))));
		return R_NilValue;
	}

	SEXP bam_alignment_SetIsSecondaryAlignment(SEXP pAlignment, SEXP bool_ok)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_SetIsSecondaryAlignment: no external pointer!");
			return R_NilValue;
		}
		if(TYPEOF(bool_ok)!=LGLSXP)
		{
			error("bamAlignment_SetIsSecondaryAlignment: no bool value!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		pAlign->SetIsSecondaryAlignment(*(LOGICAL(AS_LOGICAL(bool_ok))));
		return R_NilValue;
	}

	SEXP bam_alignment_SetIsSecondMate(SEXP pAlignment, SEXP bool_ok)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_SetIsSecondMate: no external pointer!");
			return R_NilValue;
		}
		if(TYPEOF(bool_ok)!=LGLSXP)
		{
			error("bamAlignment_SetIsSecondMate: no bool value!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		pAlign->SetIsSecondMate(*(LOGICAL(AS_LOGICAL(bool_ok))));
		return R_NilValue;
	}

	SEXP bam_alignment_SetIsUnmapped(SEXP pAlignment, SEXP bool_ok)
	{
		if(TYPEOF(pAlignment)!=EXTPTRSXP)
		{
			error("bamAlignment_SetIsUnmapped: no external pointer!");
			return R_NilValue;
		}
		if(TYPEOF(bool_ok)!=LGLSXP)
		{
			error("bamAlignment_SetIsUnmapped: no bool value!");
			return R_NilValue;
		}
		BamAlignment *pAlign=static_cast<BamAlignment*>(R_ExternalPtrAddr(pAlignment));
		pAlign->SetIsUnmapped(*(LOGICAL(AS_LOGICAL(bool_ok))));
		return R_NilValue;
	}
}

