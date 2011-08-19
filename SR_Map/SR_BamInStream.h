/*
 * =====================================================================================
 *
 *       Filename:  SR_BamInStream.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/18/2011 05:39:37 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_BAMINSTREAM_H
#define  SR_BAMINSTREAM_H

#include <stdlib.h>

#include "bam.h"
#include "SR_Types.h"


//===============================
// Type and constant definition
//===============================

// structure holding related bam input variables
typedef struct SR_BamInStreamPrvt SR_BamInStream;


//===============================
// Constructors and Destructors
//===============================

SR_BamInStream* SR_BamInStreamAlloc(const char* bamFilename);

void SR_BamInStreamFree(SR_BamInStream* pBamInStream);


//==================
// Inline functions
//==================

//===============================================================
// function:
//      see if we could jump in the bam file
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//=============================================================== 
#define SR_BamCanJump(pBamInStream) ((pBamInStream)->pBamIndex != NULL)


//======================
// Interface functions
//======================

//===============================================================
// function:
//      jump to a certain chromosome in a bam file
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. refID : the reference ID we want to jump to
// 
// return:
//      if jumping succeeds, return SR_OK; if not, return SR_ERR
//=============================================================== 
SR_Status SR_BamInStreamJump(SR_BamInStream* pBamInStream, int32_t refID);

//================================================================
// function:
//      load a unique-orphan pair from the bam file
//
// args:
//      1. ppAnchor: a pointer of pointer to the anchor mate
//      2. ppOrphan: a pointer of pointer to the orphan mate
//      2. pBamInStream : a pointer to an bam instream structure
//
// return:
//      if we get a unique-orphan pair, return SR_OK; if we reach
//      the end of file, return SR_EOF; if we finish the current
//      chromosome, return SR_OUT_OF_RANGE; else, return SR_ERR
//================================================================
SR_Status SR_BamInStreamGetPair(bam1_t** ppAnchor, bam1_t** ppOrphan, SR_BamInStream* pBamReadAux);


#endif  /*SR_BAMINSTREAM_H*/
