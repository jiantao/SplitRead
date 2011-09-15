/*
 * =====================================================================================
 *
 *       Filename:  SR_BamInStream.h
 * *    Description:  *
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
#include <gsl/gsl_histogram.h>

#include "bam.h"
#include "SR_Types.h"
#include "SR_BamMemPool.h"


//===============================
// Type and constant definition
//===============================

// structure holding related bam input variables
typedef struct SR_BamInStreamPrvt SR_BamInStream;


// stucture holding header information
typedef struct SR_BamHeader
{
    bam_header_t* pOrigHeader;

    const char** pMD5s;

}SR_BamHeader;

typedef struct SR_FragLenDstrb
{
    char** pReadGrpNames;

    void* pReadGrpHash;

    gsl_histogram** pFragLenHists;

    uint8_t drctField;

    uint32_t size;

    uint32_t capacity;

}SR_FragLenDstrb;

//===============================
// Constructors and Destructors
//===============================

SR_BamInStream* SR_BamInStreamAlloc(const char* bamFilename,        // name of input bam file
        
                                    uint32_t binLen,                // search range of a pair
                                    
                                    unsigned int numThreads,        // number of threads
                                     
                                    unsigned int buffCapacity,      // the number of alignments can be stored in each chunk of the memory pool
                                    
                                    unsigned int reportSize,        // number of unique-orphan pairs should be cached before report
                                    
                                    double scTolerance,             // soft clipping tolerance.
                                    
                                    uint8_t drctField);             // proper pair orientation flags (0 for disabling)

void SR_BamInStreamFree(SR_BamInStream* pBamInStream);

SR_BamHeader* SR_BamHeaderAlloc(void);

void SR_BamHeaderFree(SR_BamHeader* pBamHeader);

SR_FragLenDstrb* SR_FragLenDstrbAlloc(uint8_t drctField, uint32_t capacity);

void SR_FragLenDstrbFree(SR_FragLenDstrb* pDstrb);


//======================
// Inline functions
//======================

//===============================================================
// function:
//      get the number of references stored in the bam header
//
// args:
//      1. pBamHeader: a pointer to the header structure
// 
// return:
//      number of references (chromosomes)
//=============================================================== 
inline int32_t SR_BamHeaderGetRefNum(const SR_BamHeader* pBamHeader)
{
    return (pBamHeader->pOrigHeader->n_targets);
}

//===============================================================
// function:
//      get the reference ID to reference name dictionary
//
// args:
//      1. pBamHeader: a pointer to the header structure
// 
// return:
//      the dictionary of reference ID to reference name
//=============================================================== 
inline const char** SR_BamHeaderGetRefNames(const SR_BamHeader* pBamHeader)
{
    return (const char**) pBamHeader->pOrigHeader->target_name;
}

//===============================================================
// function:
//      get the array of reference length
//
// args:
//      1. pBamHeader: a pointer to the header structure
// 
// return:
//      an array contains the length of each chromosome
//=============================================================== 
inline const uint32_t* SR_BamHeaderGetRefLens(const SR_BamHeader* pBamHeader)
{
    return pBamHeader->pOrigHeader->target_len;
}

inline SR_SingleOrnt SR_BamGetOrnt(bam1_t* pAlgn)
{
    if ((pAlgn->core.flag & BAM_FREAD1) != 0)
    {
        if ((pAlgn->core.flag & BAM_FREVERSE) == 0)
            return SR_1F;
        else
            return SR_1R;
    }
    else
    {
        if ((pAlgn->core.flag & BAM_FREVERSE) == 0)
            return SR_2F;
        else
            return SR_2R;
    }
}


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
//      read the header of a bam file and load necessary
//      information from the header text
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
// 
// return:
//      if reading succeeds, return a pointer to a header
//      structure; if error is found, return NULL
//
// discussion: 
//      The file position indicator must be placed at the 
//      beginning of the file. Upon success, the position 
//      indicator will be set at the start of the first alignment
//================================================================ 
SR_BamHeader* SR_BamInStreamLoadHeader(SR_BamInStream* pBamInStream);

SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader);

void SR_FragLenDstrbSetHist(SR_FragLenDstrb* pDstrb, size_t numBins);

//================================================================
// function:
//      read an alignment from the bam file
//
// args:
//      1. pAlignment: a pointer to an alignment
//      2. pBamInStream: a pointer to an bam instream structure
// 
// return:
//      if reading succeeds, return SR_OK; if reach the end of
//      file, return SR_EOF; if get error, return SR_ERR
//
// discussion: 
//      The file position indicator must be placed at the 
//      beginning of an alignment. Upon success, the position 
//      indicator will be set at the start of the next alignment
//================================================================ 
SR_Status SR_BamInStreamRead(bam1_t* pAlignment, SR_BamInStream* pBamInStream);

//================================================================
// function:
//      get the current reference ID
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
// 
// return:
//      reference ID
//================================================================ 
int32_t SR_BamInStreamGetRefID(const SR_BamInStream* pBamInStream);

//==================================================================
// function:
//      load a certain number of unique-orphan pairs from the 
//      bam file into the buffer for a give thread
//
// args:
//      1. pBamInStream : a pointer to an bam instream structure
//      2. threadID: the ID of a thread
//      3. pDstb: a pointer to the fragment length distribution
//                object
//
// return:
//      if we get enough unique-orphan pair, return SR_OK; 
//      if we reach the end of file, return SR_EOF; if we finish 
//      the current chromosome, return SR_OUT_OF_RANGE; 
//      else, return SR_ERR
//==================================================================
SR_Status SR_BamInStreamLoadPairs(SR_BamInStream* pBamInStream, unsigned int threadID, SR_FragLenDstrb* pDstrb);

//================================================================
// function:
//      get a iterator to a certain buffer of a thread
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. threadID: ID of a thread
// 
// return:
//      iterator to the buffer of a thread
//================================================================ 
SR_BamListIter SR_BamInStreamGetIter(SR_BamInStream* pBamInStream, unsigned int threadID);

SR_BamList* SR_BamInStreamGetBuff(SR_BamInStream* pBamInStream, unsigned int threadID);

//================================================================
// function:
//      clear the alignment buffer that associates with a
//      certain thread
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. threadID: ID of a thread
//================================================================ 
void SR_BamInStreamClearBuff(SR_BamInStream* pBamInStream, unsigned int threadID);

//================================================================
// function:
//      get the size of the memory pool in a bam in stream object
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
// 
// return:
//      number of memory buffers in the memory pool
//================================================================ 
unsigned int SR_BamInStreamGetPoolSize(SR_BamInStream* pBamInStream);

//================================================================
// function:
//      decrease the size of memory pool inside the bam in stream
//      object to save memory
//
// args:
//      1. pBamInStream : a pointer to an bam instream structure
//      2. newSize: the new size of the memory pool that user 
//                  want to set (which can be only smaller than
//                  current size otherwise nothing will be done)
//
// return:
//      the actual size of the memory pool after shrinking
//
// discussion:
//      this function should be called after the processing of 
//      a certain chromosome. The memory allocated for the 
//      return lists should be freed before you can call this
//      function. the desired size may not be achieved since
//      there may not be enough empty buffer chunks. check the
//      return value for the actual size of the memory pool
//================================================================
unsigned int SR_BamInStreamShrinkPool(SR_BamInStream* pBamInStream, unsigned int newSize);

#endif  /*SR_BAMINSTREAM_H*/
