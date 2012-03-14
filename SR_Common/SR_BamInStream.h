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


#include "bam.h"
#include "SR_Types.h"
#include "SR_BamHeader.h"
#include "SR_BamMemPool.h"

//===============================
// Type and constant definition
//===============================

typedef SR_StreamCode (*SR_BamFilter) (const bam1_t* pAlignment, void* pFilterData, int32_t currRefID, int32_t currBinPos);

// control parameters for bam in stream
typedef enum SR_StreamControlFlag
{
    SR_NO_SPECIAL_CONTROL = 0,

    SR_USE_BAM_INDEX      = 1,   // bam in stream will open the bam index file and load it into memory

    SR_READ_PAIR_MODE     = 2    // bam in stream will find read pairs genomically

}SR_StreamControlFlag;

typedef struct SR_StreamMode
{
    SR_BamFilter filterFunc;             // a filter function used to skip those uninterested reads

    void* filterData;                    // parameters for the filter function

    SR_StreamControlFlag controlFlag;    // flag used to control the bam in stream

}SR_StreamMode;

// alignment type
typedef enum
{
    SR_UNIQUE_ORPHAN = 0,      // one mate is uniquely aligned the other mate is unaligned

    SR_UNIQUE_SOFT = 1,        // one mate is uniquely aligned the other mate has may soft-clipping bases

    SR_UNIQUE_MULTIPLE = 2,    // one mate is uniquely aligned the other mate is mutiply aligned

    SR_UNIQUE_NORMAL = 3,      // both mates are uniquely aligned

    SR_OTHER_ALGN_TYPE = 4     // other alignment types

}SR_AlgnType;

// iterator used to retrieve the alignments and alignment types
typedef struct SR_BamInStreamIter
{
    SR_BamNode* pBamNode;    // a pointer to the bam node structure

    SR_AlgnType* pAlgnType;  // a pointer to the alignment type

}SR_BamInStreamIter;

// private data structure that holds all bam-input-related information
typedef struct SR_BamInStream
{
    bamFile fpBamInput;                        // file pointer to a input bam file

    bam_index_t* pBamIndex;                    // file pointer to a input bam index file

    SR_BamFilter filterFunc;                   // customized filter function 

    void* filterData;                          // data used by the filter function

    SR_BamMemPool* pMemPool;                   // memory pool used to allocate and recycle the bam alignments

    void* pNameHashes[2];                      // two hashes used to get a pair of alignments

    SR_BamList* pRetLists;                     // when we find any unique-orphan pairs we push them into these lists, each thread has its own list

    SR_AlgnType* pAlgnTypes;                   // store the alignment types of read pairs(unique-orphan, unique-multiple, unique-softclipping, unique-unique...)

    SR_BamNode* pNewNode;                      // the just read-in bam alignment

    SR_BamList pAlgnLists[2];                  // lists used to store those incoming alignments

    unsigned int numThreads;                   // number of threads will be used

    unsigned int reportSize;                   // number of alignments should be loaded before report

    int32_t currRefID;                         // the reference ID of the current read-in alignment

    int32_t currBinPos;                        // the start position of current bin (0-based)

    uint32_t binLen;                           // the length of the bin

    SR_StreamControlFlag controlFlag; 

}SR_BamInStream;

typedef struct SR_SplitBamInStream
{
    bamFile fpBamInput;

    int32_t currRefID;

}SR_SplitBamInStream;


//===============================
// Constructors and Destructors
//===============================

SR_BamInStream* SR_BamInStreamAlloc(uint32_t binLen,                       // search range of a pair
                                    
                                    unsigned int numThreads,               // number of threads
                                     
                                    unsigned int buffCapacity,             // the number of alignments can be stored in each chunk of the memory pool
                                    
                                    unsigned int reportSize,               // number of alignments should be cached before report
                                    
                                    SR_StreamMode* pStreamMode);           // a structure used to control the features of the stream


void SR_BamInStreamFree(SR_BamInStream* pBamInStream);

SR_SplitBamInStream* SR_SplitBamInStreamAlloc(void);

void SR_SplitBamInStreamFree(SR_SplitBamInStream* pSplitBamInStream);


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

//===============================================================
// function:
//      open a bam file
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. bamFileName: the name of the bam file
// 
// return:
//      if open succeeds, return SR_OK; if not, return SR_ERR
//=============================================================== 
SR_Status SR_BamInStreamOpen(SR_BamInStream* pBamInStream, const char* bamFileName);

//===============================================================
// function:
//      clear the bam instream object(return list, 
//      alignment list and name hash)
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//=============================================================== 
void SR_BamInStreamClear(SR_BamInStream* pBamInStream);

//===============================================================
// function:
//      close the current bam files and clear the bam instream
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//=============================================================== 
void SR_BamInStreamClose(SR_BamInStream* pBamInStream);

//===============================================================
// function:
//      tell the current position in the bam file
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//
// return:
//      the current position in the bam file
//=============================================================== 
#define SR_BamInStreamTell(pBamInStream) bam_tell((pBamInStream)->fpBamInput)

//===============================================================
// function:
//      seek to a specific position in the bam file
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. pos: file position
//      3. where: from where we should seek 
//                (only SEEK_SET is allowed)
//=============================================================== 
#define SR_BamInStreamSeek(pBamInStream, pos, where) bam_seek((pBamInStream)->fpBamInput, pos, where)

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

//===============================================================
// function:
//      set the mode of the bam in stream
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. filterFunc: filter function for the bam in stream
//      3. filterData: filter data passed to filter function
//      4. controlFlag: flags used to control the stream
// 
//=============================================================== 
static inline void SR_SetStreamMode(SR_StreamMode* pStreamMode, SR_BamFilter filterFunc, void* filterData, SR_StreamControlFlag controlFlag)
{
    pStreamMode->filterFunc = filterFunc;
    pStreamMode->filterData = filterData;
    pStreamMode->controlFlag = controlFlag;
}

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
static inline SR_Status SR_BamInStreamRead(bam1_t* pAlignment, SR_BamInStream* pBamInStream)
{
    int ret = bam_read1(pBamInStream->fpBamInput, pAlignment);

    if (ret > 0)
        return SR_OK;
    else if (ret == -1)
        return SR_EOF;
    else
        return SR_ERR;
}

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
#define SR_BamInStreamGetRefID(pBamInStream) ((pBamInStream)->currRefID)

//==================================================================
// function:
//      load a pair of bam alignments
//
// args:
//      1. ppAlgnOne: a pointer to the pointer of an alignment
//      2. ppAlgnTwo: a pointer to the pointer of an alignment
//      3. pBamInStream : a pointer to an bam instream structure
//
// return:
//      if we get enough unique-orphan pair, return SR_OK; 
//      if we reach the end of file, return SR_EOF; if we finish 
//      the current chromosome, return SR_OUT_OF_RANGE; 
//      else, return SR_ERR
//==================================================================
SR_Status SR_BamInStreamLoadPair(SR_BamNode** ppAlgnOne, SR_BamNode** ppAlgnTwo, SR_BamInStream* pBamInStream);

//==================================================================
// function:
//      get the alignment type from a read pair
//
// args:
//      ppAnchor: a pointer to a pointer of bam node with anchor
//                alignment
//      ppOrphan: a pointer to a pointer of bam node with orphan
//                alignment
//      scTolerance: soft clipping tolerance
//
// return:
//      the alignment type of the read pair
//==================================================================
SR_AlgnType SR_GetAlignmentType(SR_BamNode** ppAlgnOne, SR_BamNode** ppAlgnTwo, double scTolerance, double maxMismatchRate, unsigned char minMQ);

//==================================================================
// function:
//      load a certain number of unique orphan pairs into a buffer
//      associated with a thread
//
// args:
//      1. pBamInStream: a pointer to a bam in stream object
//      2. threadID: the ID of the thread
//      3. scTolerance: soft clipping tolerance
//
// return:
//      status of the bam in stream. if we reach the end of a
//      chromosome, return SR_OUT_OF_RANGE; if we reach the end of
//      the file, return SR_EOF; if an error happens, return
//      SR_ERR; else return SR_OK
//==================================================================
SR_Status SR_LoadAlgnPairs(SR_BamInStream* pBamInStream, unsigned int threadID, double scTolerance, double maxMismatchRate, unsigned char minMQ);

//================================================================
// function:
//      get the size of the memory pool in the bam in stream
//      structure
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
// 
// return:
//      the size of memory pool in the bam in stream object
//================================================================ 
#define SR_BamInStreamGetPoolSize(pBamInStream) ((pBamInStream)->pMemPool->numBuffs)

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
#define SR_BamInStreamSetIter(pIter, pBamInStream, threadID)                                           \
    do                                                                                                 \
    {                                                                                                  \
        (pIter)->pBamNode = SR_BamListGetIter((pBamInStream)->pRetLists + (threadID));                 \
        (pIter)->pAlgnType = (pBamInStream)->pAlgnTypes + (pBamInStream)->reportSize * (threadID);     \
                                                                                                       \
    }while(0)

#define SR_BamInStreamSetAlgnType(pBamInStream, threadID, algnType)                                                            \
    do                                                                                                                         \
    {                                                                                                                          \
        unsigned int index = (pBamInStream)->reportSize * (threadID) + (pBamInStream)->pRetLists[(threadID)].numNode / 2 -1;   \
        (pBamInStream)->pAlgnTypes[index] = (algnType);                                                                        \
                                                                                                                               \
    }while(0)

//================================================================
// function:
//      push the qualified alignment into a given thread buffer
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. pAlignment: a pointer to a qualified alignment
//      2. threadID: ID of a thread
// 
// return:
//      status of the thread buffer. if the thread buffer is full
//      then return SR_FULL else return SR_OK
//================================================================ 
static inline SR_Status SR_BamInStreamPush(SR_BamInStream* pBamInStream, SR_BamNode* pAlignment, unsigned int threadID)
{
    SR_BamListPushBack(pBamInStream->pRetLists + threadID, pAlignment);

    if (pBamInStream->pRetLists[threadID].numNode == pBamInStream->reportSize)
        return SR_FULL;

    return SR_OK;
}

//================================================================
// function:
//      recycle an unwanted bam node
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. pBamNode: a pointer to a bam node 
//================================================================ 
#define SR_BamInStreamRecycle(pBamInStream, pBamNode) SR_BamListPushHead(&((pBamInStream)->pMemPool->avlNodeList), (pBamNode))

//================================================================
// function:
//      clear the return list  associated with a certain thread
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. threadID: the id of the thread that needs to be clear
//================================================================ 
#define SR_BamInStreamClearRetList(pBamInStream, threadID) SR_BamListReset((pBamInStream)->pRetLists + (threadID), (pBamInStream)->pMemPool)

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


void SR_SplitBamInStreamOpen(SR_SplitBamInStream* pSplitBamInStream, const char* fileName);

void SR_SplitBamInStreamClose(SR_SplitBamInStream* pSplitBamInStream);

SR_Status SR_SplitBamInStreamRead(bam1_t (*pAlgns)[3], SR_SplitBamInStream* pSplitBamInStream);

#endif  /*SR_BAMINSTREAM_H*/
