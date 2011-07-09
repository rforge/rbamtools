/*
 * rbamtools.h
 */

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>

#include <cstdio>
#include <string>
#include <vector>

#include "bamtools/BamReader.h"
#include "bamtools/BamWriter.h"

using namespace std;
using namespace BamTools;

extern "C"{

///////////////////////////////////////////////////////////////////////////////////////////////
// BamWriter
///////////////////////////////////////////////////////////////////////////////////////////////
static void finalize_bam_writer(SEXP ptr);
SEXP get_bam_writer();
SEXP bam_writer_open_reader(SEXP pWriter, SEXP pReader,SEXP pFilename);
vector<RefData> bam_writer_dfToVector(SEXP pRefData);
SEXP bam_writer_open(SEXP pWriter,SEXP pFilename,SEXP pSamHeader,SEXP pRefSeqs);
SEXP bam_writer_SaveAlignment(SEXP pWriter, SEXP pAlignment);
SEXP bam_writer_close(SEXP pWriter);


///////////////////////////////////////////////////////////////////////////////////////////////
// BamReader
///////////////////////////////////////////////////////////////////////////////////////////////
static void finalize_bam_reader(SEXP ptr);
SEXP get_bam_reader();
SEXP bam_reader_open(SEXP pReader,SEXP filename,SEXP index_file_name);
SEXP bam_reader_close(SEXP pReader);
SEXP bam_reader_GetHeaderText(SEXP pReader);
SEXP bam_reader_GetReferenceCount(SEXP pReader);
SEXP bam_reader_GetReferenceData(SEXP pReader);
SEXP bam_reader_GetReferenceID(SEXP pReader, SEXP refName);
SEXP bam_reader_Jump(SEXP pReader,SEXP refID,SEXP position);
SEXP bam_reader_Rewind(SEXP pReader);
SEXP bam_reader_CreateIndex(SEXP pReader);
SEXP bam_reader_GetNextAlignment(SEXP pReader);

///////////////////////////////////////////////////////////////////////////////////////////////
// BamAlignment
///////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_BamAlignment(SEXP pAlignment);
SEXP bam_alignment_getName(SEXP pAlignment);
SEXP bam_alignment_getRefID(SEXP pAlignment);
SEXP bam_alignment_getPosition(SEXP pAlignment);
SEXP bam_alignment_getCigar_df(SEXP pAlignment);
SEXP bam_alignment_getMateRefID(SEXP pAlignment);
SEXP bam_alignment_getMatePosition(SEXP pAlignment);
SEXP bam_alignment_getInsertSize(SEXP pAlignment);
SEXP bam_alignment_getMapQuality(SEXP pAlignment);
SEXP bam_alignment_getQueryBases(SEXP pAlignment);
SEXP bam_alignment_getAlignedBases(SEXP pAlignment);
SEXP bam_alignment_getQualities(SEXP pAlignment);

/////////////////////////////////////////////////
// BamAux.h queries against alignment flag

// Reading queries
SEXP bam_alignment_IsDuplicate(SEXP pAlignment);
SEXP bam_alignment_IsFailedQC(SEXP pAlignment);
SEXP bam_alignment_IsFirstMate(SEXP pAlignment);
SEXP bam_alignment_IsMapped(SEXP pAlignment);
SEXP bam_alignment_IsMateMapped(SEXP pAlignment);
SEXP bam_alignment_IsMateReverseStrand(SEXP pAlignment);
SEXP bam_alignment_IsPaired(SEXP pAlignment);
SEXP bam_alignment_IsPrimaryAlignment(SEXP pAlignment);
SEXP bam_alignment_IsProperPair(SEXP pAlignment);
SEXP bam_alignment_IsReverseStrand(SEXP pAlignment);
SEXP bam_alignment_IsSecondMate(SEXP pAlignment);

// Writing queries
SEXP bam_alignment_SetIsDuplicate(SEXP pAlignment, SEXP bool_ok);
SEXP bam_alignment_SetIsFailedQC(SEXP pAlignment, SEXP bool_ok);
SEXP bam_alignment_SetIsFirstMate(SEXP pAlignment, SEXP bool_ok);
SEXP bam_alignment_SetIsMateUnmapped(SEXP pAlignment, SEXP bool_ok);
SEXP bam_alignment_SetIsMateReverseStrand(SEXP pAlignment, SEXP bool_ok);
SEXP bam_alignment_SetIsPaired(SEXP pAlignment, SEXP bool_ok);
SEXP bam_alignment_SetIsProperPair(SEXP pAlignment, SEXP bool_ok);
SEXP bam_alignment_SetIsReverseStrand(SEXP pAlignment, SEXP bool_ok);
SEXP bam_alignment_SetIsSecondaryAlignment(SEXP pAlignment, SEXP bool_ok);
SEXP bam_alignment_SetIsSecondMate(SEXP pAlignment, SEXP bool_ok);
SEXP bam_alignment_SetIsUnmapped(SEXP pAlignment, SEXP bool_ok);

}
