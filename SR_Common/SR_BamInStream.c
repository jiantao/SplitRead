/*
 * =====================================================================================
 *
 *       Filename:  SR_BamInStream.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/18/2011 05:41:34 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <assert.h>

#include "khash.h"
#include "SR_Error.h"
#include "SR_Utilities.h"
#include "SR_BamInStream.h"


//===============================
// Type and constant definition
//===============================

// value indicate that we did not load any read
// into the bam array
#define NO_QUERY_YET (-2)

// default capacity of a bam array
#define DEFAULT_BAM_ARRAY_CAP 200

#define SR_MAX_BIN_LEN 500000000

// a mask used to filter out those unwanted reads for split alignments
// it includes proper paired reads, secondar reads, qc-failed reads and duplicated reads
#define SR_BAM_FMASK (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)

// bin number of buffers
#define PREV_BIN 0
#define CURR_BIN 1

// default soft clipping tolerance
#define DEFAULT_SC_TOLERANCE 0.2

// default mismatch tolerance
#define DEFAULT_MAX_MISMATCH_RATE 0.1


// alignment status
enum AlignmentStatus
{
    NONE_GOOD   = -1,       // neither a good anchor nor a good orphan candidate

    GOOD_ANCHOR    = 0,     // a good anchor candidate

    GOOD_ORPHAN    = 1,     // a good orphan candidate

    GOOD_SOFT      = 2,     // a good soft clipping candidate

    GOOD_MULTIPLE  = 3,     // a good multiple aligned candidate
};

// initialize a string-to-bam hash table used to retrieve a pair of read
KHASH_MAP_INIT_STR(queryName, SR_BamNode*);

KHASH_SET_INIT_INT64(buffAddress);


//===================
// Static functions
//===================

static inline int SR_BamInStreamLoadNext(SR_BamInStream* pBamInStream)
{
    // for the bam alignment array, if we need to expand its space
    // we have to initialize those newly created bam alignment 
    // and update the query name hash since the address of those
    // bam alignments are changed after expanding
    pBamInStream->pNewNode = SR_BamNodeAlloc(pBamInStream->pMemPool);
    if (pBamInStream->pNewNode == NULL)
        SR_ErrQuit("ERROR: Too many unpaired reads are stored in the memory. Please use smaller bin size or disable searching pair genomically.\n");

    int ret = bam_read1(pBamInStream->fpBamInput, &(pBamInStream->pNewNode->alignment));

    return ret;
}

static void SR_BamInStreamReset(SR_BamInStream* pBamInStream)
{
    pBamInStream->pNewNode = NULL;

    pBamInStream->currBinPos = NO_QUERY_YET;
    pBamInStream->currRefID = NO_QUERY_YET;

    kh_clear(queryName, pBamInStream->pNameHashes[PREV_BIN]);
    kh_clear(queryName, pBamInStream->pNameHashes[CURR_BIN]);

    SR_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);
    SR_BamListReset(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pMemPool);
}

static double SR_GetMismatchRate(const bam1_t* pAlignment)
{
    uint32_t* cigar = bam1_cigar(pAlignment);
    unsigned int numMismatch = 0;
    for (unsigned i = 0; i != pAlignment->core.n_cigar; ++i)
    {
        int type = (cigar[i] & BAM_CIGAR_MASK);
        if (type == BAM_CINS || type == BAM_CDEL || type == BAM_CSOFT_CLIP || type == BAM_CMISMATCH)
        {
            numMismatch += (cigar[i] >> BAM_CIGAR_SHIFT);
        }
    }

    return ((double) numMismatch / pAlignment->core.l_qseq);
}

static int SR_CheckAlignment(const bam1_t* pAlignment, double scTolerance, double maxMismatchRate, unsigned char minMQ)
{
    if ((pAlignment->core.flag & BAM_FUNMAP) != 0)
        return GOOD_ORPHAN;

    unsigned int scLimit = scTolerance * pAlignment->core.l_qseq;
    uint32_t* cigar = bam1_cigar(pAlignment);

    SR_Bool isHeadSC = FALSE;
    SR_Bool isTailSC = FALSE;

    if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP 
        && (cigar[0] >> BAM_CIGAR_SHIFT) >= scLimit)
    {
        isHeadSC = TRUE;
    }

    unsigned int lastIndex = pAlignment->core.n_cigar - 1;
    if ((cigar[lastIndex] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP 
        && (cigar[lastIndex] >> BAM_CIGAR_SHIFT) >= scLimit)
    {
        isTailSC = TRUE;
    }

    if (!isHeadSC && !isTailSC) 
    {
        if (pAlignment->core.qual >= minMQ)
            return GOOD_ANCHOR;
        else if (maxMismatchRate > 0.0)
        {
            double mismatchRate = SR_GetMismatchRate(pAlignment);
            if (mismatchRate <= maxMismatchRate)
                return GOOD_MULTIPLE;
        }

        return NONE_GOOD;
    }
    else if (isHeadSC && isTailSC)
        return NONE_GOOD;
    else
        return GOOD_SOFT;
}

//===============================
// Constructors and Destructors
//===============================

SR_BamInStream* SR_BamInStreamAlloc( uint32_t binLen, unsigned int numThreads, unsigned int buffCapacity, 
                                    unsigned int reportSize, SR_StreamMode* pStreamMode)
{
    SR_BamInStream* pBamInStream = (SR_BamInStream*) calloc(1, sizeof(SR_BamInStream));
    if (pBamInStream == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a bam input stream object.");

    pBamInStream->fpBamInput = NULL;
    pBamInStream->pBamIndex = NULL;

    pBamInStream->filterFunc = pStreamMode->filterFunc;
    pBamInStream->filterData = pStreamMode->filterData;
    pBamInStream->controlFlag = pStreamMode->controlFlag;

    pBamInStream->numThreads = numThreads;
    pBamInStream->reportSize = reportSize;

    pBamInStream->currRefID = NO_QUERY_YET;
    pBamInStream->currBinPos = NO_QUERY_YET;
    pBamInStream->binLen = binLen;
    pBamInStream->pNewNode = NULL;

    if (numThreads > 0)
    {
        pBamInStream->pRetLists = (SR_BamList*) calloc(numThreads, sizeof(SR_BamList));
        if (pBamInStream->pRetLists == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of retrun alignment lists in the bam input stream object.\n");

        pBamInStream->pAlgnTypes = (SR_AlgnType*) malloc(numThreads * reportSize * sizeof(SR_AlgnType));
        if (pBamInStream->pAlgnTypes == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of pair alignment type in the bam input stream object.\n");
    }
    else
    {
        pBamInStream->pRetLists = NULL;
        pBamInStream->pAlgnTypes = NULL;
        pBamInStream->reportSize = 0;
    }

    pBamInStream->pNameHashes[PREV_BIN] = kh_init(queryName);
    if (pBamInStream->reportSize > 0)
        kh_resize(queryName, pBamInStream->pNameHashes[PREV_BIN], reportSize);

    pBamInStream->pNameHashes[CURR_BIN] = kh_init(queryName);
    if (pBamInStream->reportSize > 0)
        kh_resize(queryName, pBamInStream->pNameHashes[CURR_BIN], reportSize);

    pBamInStream->pMemPool = SR_BamMemPoolAlloc(buffCapacity);

    return pBamInStream;
}

void SR_BamInStreamFree(SR_BamInStream* pBamInStream)
{
    if (pBamInStream != NULL)
    {
        kh_destroy(queryName, pBamInStream->pNameHashes[PREV_BIN]);
        kh_destroy(queryName, pBamInStream->pNameHashes[CURR_BIN]);

        free(pBamInStream->pRetLists);
        free(pBamInStream->pAlgnTypes);
        SR_BamMemPoolFree(pBamInStream->pMemPool);

        free(pBamInStream);
    }
}


//======================
// Interface functions
//======================

SR_Status SR_BamInStreamOpen(SR_BamInStream* pBamInStream, const char* bamFilename)
{
    pBamInStream->fpBamInput = bam_open(bamFilename, "r");
    if (pBamInStream->fpBamInput == NULL)
    {
        SR_ErrMsg("ERROR: Cannot open bam file: %s.\n", bamFilename);
        return SR_ERR;
    }

    if ((pBamInStream->controlFlag & SR_USE_BAM_INDEX) != 0)
    {
        pBamInStream->pBamIndex = bam_index_load(bamFilename);
        if (pBamInStream->pBamIndex == NULL)
        {
            SR_ErrMsg("WARNING: Cannot open bam index file for: %s\n", bamFilename);
            return SR_ERR;
        }
    }

    return SR_OK;
}

void SR_BamInStreamClear(SR_BamInStream* pBamInStream)
{
    pBamInStream->pNewNode = NULL;
    pBamInStream->currRefID = NO_QUERY_YET;
    pBamInStream->currBinPos = NO_QUERY_YET;

    for (unsigned int i = 0; i != pBamInStream->numThreads; ++i)
        SR_BamListReset(pBamInStream->pRetLists + i, pBamInStream->pMemPool);

    for (unsigned int i = 0; i != 2; ++i)
        SR_BamListReset(pBamInStream->pAlgnLists + i, pBamInStream->pMemPool);

    kh_clear(queryName, pBamInStream->pNameHashes[PREV_BIN]);
    kh_clear(queryName, pBamInStream->pNameHashes[CURR_BIN]);
}

void SR_BamInStreamClose(SR_BamInStream* pBamInStream)
{
    bam_close(pBamInStream->fpBamInput);
    pBamInStream->fpBamInput = NULL;

    if (pBamInStream->pBamIndex != NULL)
    {
        bam_index_destroy(pBamInStream->pBamIndex);
        pBamInStream->pBamIndex = NULL;
    }

    SR_BamInStreamClear(pBamInStream);
}

// jump to a certain chromosome in a bam file
SR_Status SR_BamInStreamJump(SR_BamInStream* pBamInStream, int32_t refID)
{
    // if we do not have the index file return error
    if (pBamInStream->pBamIndex == NULL)
        return SR_ERR;

    // clear the bam array before jump
    SR_BamInStreamReset(pBamInStream);

    // jump and read the first alignment in the given chromosome
    int ret;
    bam_iter_t pBamIter = bam_iter_query(pBamInStream->pBamIndex, refID, 0, INT_MAX);

    pBamInStream->pNewNode = SR_BamNodeAlloc(pBamInStream->pMemPool);
    ret = bam_iter_read(pBamInStream->fpBamInput, pBamIter, &(pBamInStream->pNewNode->alignment));
    bam_iter_destroy(pBamIter);

    khash_t(queryName)* pNameHashCurr = NULL;

    // see if we jump to the desired chromosome
    if (ret > 0 && pBamInStream->pNewNode->alignment.core.tid == refID)
    {
        // exclude those reads who are non-paired-end, qc-fail, duplicate-marked, proper-paired, 
        // both aligned, secondary-alignment and no-name-specified.
        if ((pBamInStream->pNewNode->alignment.core.flag & SR_BAM_FMASK) != 0
            || strcmp(bam1_qname(&(pBamInStream->pNewNode->alignment)), "*") == 0)
        {
            SR_BamNodeFree(pBamInStream->pNewNode, pBamInStream->pMemPool);
            pBamInStream->pNewNode = NULL;
            pBamInStream->currBinPos = NO_QUERY_YET;
        }
        else
        {
            SR_BamListPushHead(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pNewNode);

            int khRet = 0;
            khiter_t khIter = kh_put(queryName, pBamInStream->pNameHashes[CURR_BIN], bam1_qname(&(pBamInStream->pNewNode->alignment)), &khRet);

            if (khRet != 0)
            {
                pNameHashCurr = pBamInStream->pNameHashes[CURR_BIN];
                kh_value(pNameHashCurr, khIter) = pBamInStream->pNewNode;
            }
            else
                return SR_ERR;

            pBamInStream->currBinPos = pBamInStream->pNewNode->alignment.core.pos; 
            pBamInStream->pNewNode = NULL;
        }

        pBamInStream->currRefID = refID;
        return SR_OK;
    }
    else if (ret == -1)
    {
        return SR_OUT_OF_RANGE;
    }
    else
    {
        return SR_ERR;
    }
}

// read the header of a bam file
SR_BamHeader* SR_BamInStreamLoadHeader(SR_BamInStream* pBamInStream)
{
    bam_header_t* pOrigHeader = bam_header_read(pBamInStream->fpBamInput);
    if (pOrigHeader == NULL)
        return NULL;

    SR_BamHeader* pBamHeader = SR_BamHeaderAlloc();

    pBamHeader->pOrigHeader = pOrigHeader;


    pBamHeader->pMD5s = (const char**) calloc(pOrigHeader->n_targets, sizeof(char*));
    if (pBamHeader->pMD5s == NULL)
        SR_ErrQuit("ERROR: Not enough memory for md5 string");

    unsigned int numMD5 = 0;
    for (const char* md5Pos = pOrigHeader->text; numMD5 <= pOrigHeader->n_targets && (md5Pos = strstr(md5Pos, "M5:")) != NULL; ++numMD5, ++md5Pos)
    {
        pBamHeader->pMD5s[numMD5] = md5Pos + 3;
    }

    if (numMD5 != pOrigHeader->n_targets)
    {
        free(pBamHeader->pMD5s);
        pBamHeader->pMD5s = NULL;

        if (numMD5 != 0)
            SR_ErrMsg("WARNING: Number of MD5 string is not consistent with number of chromosomes.");
    }

    return pBamHeader;
}

void SR_CheckBamInStream(int* pPrevHashSize, int* pCurrHashSize, int* pPrevListSize, int* pCurrListSize, const SR_BamInStream* pBamInStream)
{
    khash_t(queryName)* pNameHashPrev = pBamInStream->pNameHashes[PREV_BIN];
    khash_t(queryName)* pNameHashCurr = pBamInStream->pNameHashes[CURR_BIN];

    *pPrevHashSize = kh_size(pNameHashPrev);
    *pCurrHashSize = kh_size(pNameHashCurr);

    *pPrevListSize = pBamInStream->pAlgnLists[PREV_BIN].numNode;
    *pCurrListSize = pBamInStream->pAlgnLists[CURR_BIN].numNode;
}

// load a unique-orphan pair from a bam file
SR_Status SR_BamInStreamLoadPair(SR_BamNode** ppUpAlgn, SR_BamNode** ppDownAlgn, SR_BamInStream* pBamInStream) 
{
    (*ppUpAlgn) = NULL;
    (*ppDownAlgn) = NULL;

    khash_t(queryName)* pNameHashPrev = pBamInStream->pNameHashes[PREV_BIN];
    khash_t(queryName)* pNameHashCurr = pBamInStream->pNameHashes[CURR_BIN];

    int ret = 1;
    while(ret > 0 && (ret = SR_BamInStreamLoadNext(pBamInStream)) > 0)
    {
        // exclude those reads who are non-paired-end, qc-fail, duplicate-marked, proper-paired, 
        // both aligned, secondary-alignment and no-name-specified.
        SR_StreamCode filterCode = pBamInStream->filterFunc(&(pBamInStream->pNewNode->alignment), pBamInStream->filterData, 
                                                            pBamInStream->currRefID, pBamInStream->currBinPos);

        if (filterCode != STREAM_KEEP)
        {
            SR_BamNodeFree(pBamInStream->pNewNode, pBamInStream->pMemPool);
            pBamInStream->pNewNode = NULL;

            if (filterCode == STREAM_PASS)
                continue;
            else
            {
                ret = SR_OK;
                break;
            }
        }

        // update the current ref ID or position if the incoming alignment has a 
        // different value. The name hash and the bam array will be reset
        if (pBamInStream->pNewNode->alignment.core.pos >= pBamInStream->currBinPos + 2 * pBamInStream->binLen
            || pBamInStream->pNewNode->alignment.core.tid != pBamInStream->currRefID
            || pBamInStream->currBinPos == NO_QUERY_YET)
        {
            if (pBamInStream->pNewNode->alignment.core.tid != pBamInStream->currRefID
                && pBamInStream->currRefID != NO_QUERY_YET)
            {
                ret = SR_OUT_OF_RANGE;
            }

            pBamInStream->currRefID  = pBamInStream->pNewNode->alignment.core.tid;
            pBamInStream->currBinPos = pBamInStream->pNewNode->alignment.core.pos;

            kh_clear(queryName, pNameHashPrev);
            kh_clear(queryName, pNameHashCurr);

            SR_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);
            SR_BamListReset(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pMemPool);

        }
        else if (pBamInStream->pNewNode->alignment.core.pos >= pBamInStream->currBinPos + pBamInStream->binLen)
        {
            pBamInStream->currBinPos += pBamInStream->binLen;

            kh_clear(queryName, pNameHashPrev);
            SR_SWAP(pNameHashPrev, pNameHashCurr, khash_t(queryName)*);

            SR_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);

            SR_SWAP(pBamInStream->pAlgnLists[PREV_BIN], pBamInStream->pAlgnLists[CURR_BIN], SR_BamList);
        }

        SR_BamListPushHead(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pNewNode);

        (*ppUpAlgn) = NULL;
        (*ppDownAlgn) = NULL;

        khiter_t khIter = 0;

        khIter = kh_get(queryName, pNameHashPrev, bam1_qname(&(pBamInStream->pNewNode->alignment)));

        if (khIter != kh_end(pNameHashPrev))
        {
            ret = SR_OK;
            (*ppUpAlgn) = kh_value(pNameHashPrev, khIter);
            (*ppDownAlgn) = pBamInStream->pNewNode;

            kh_del(queryName, pNameHashPrev, khIter);

            SR_BamListRemove(&(pBamInStream->pAlgnLists[PREV_BIN]), (*ppUpAlgn));
            SR_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), (*ppDownAlgn));
        }
        else
        {
            int khRet = 0;
            khIter = kh_put(queryName, pNameHashCurr, bam1_qname(&(pBamInStream->pNewNode->alignment)), &khRet);

            if (khRet == 0) // we found a pair of alignments 
            {
                ret = SR_OK;
                (*ppUpAlgn) = kh_value(pNameHashCurr, khIter);
                (*ppDownAlgn) = pBamInStream->pNewNode;

                kh_del(queryName, pNameHashCurr, khIter);

                SR_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), (*ppUpAlgn));
                SR_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), (*ppDownAlgn));
            }
            else // not finding corresponding mate, save the current value and move on
            {
                kh_value(pNameHashCurr, khIter) = pBamInStream->pNewNode;
            }
        }
    }

    pBamInStream->pNameHashes[PREV_BIN] = pNameHashPrev;
    pBamInStream->pNameHashes[CURR_BIN] = pNameHashCurr;

    if (ret < 0)
    {
        if (ret != SR_OUT_OF_RANGE)
            SR_BamNodeFree(pBamInStream->pNewNode, pBamInStream->pMemPool);

        if ( ret != SR_OUT_OF_RANGE && ret != SR_EOF)
            ret = SR_ERR;
    }

    pBamInStream->pNewNode = NULL;
    return ret;
}

SR_AlgnType SR_GetAlignmentType(SR_BamNode** ppAlgnOne, SR_BamNode** ppAlgnTwo, double scTolerance, double maxMismatchRate, unsigned char minMQ)
{
    int firstType = SR_CheckAlignment(&((*ppAlgnOne)->alignment), scTolerance, maxMismatchRate, minMQ);
    int secondType = SR_CheckAlignment(&((*ppAlgnTwo)->alignment), scTolerance, maxMismatchRate, minMQ);

    if (firstType != GOOD_ANCHOR && secondType == GOOD_ANCHOR)
    {
        SR_SWAP(*ppAlgnOne, *ppAlgnTwo, SR_BamNode*);
        SR_SWAP(firstType, secondType, int);
    }

    if (firstType == GOOD_ANCHOR)
    {
        switch(secondType)
        {
            case GOOD_ORPHAN:
                return SR_UNIQUE_ORPHAN;
                break;
            case GOOD_SOFT:
                return SR_UNIQUE_SOFT;
                break;
            case GOOD_MULTIPLE:
                return SR_UNIQUE_MULTIPLE;
                break;
            case GOOD_ANCHOR:
                return SR_UNIQUE_NORMAL;
                break;
            default:
                return SR_OTHER_ALGN_TYPE;
                break;
        }
    }

    return SR_OTHER_ALGN_TYPE;
}

SR_Status SR_LoadAlgnPairs(SR_BamInStream* pBamInStream, unsigned int threadID, double scTolerance, double maxMismatchRate, unsigned char minMQ)
{
    SR_BamNode* pAlgnOne = NULL;
    SR_BamNode* pAlgnTwo = NULL;

    SR_Status readerStatus = SR_OK;
    SR_Status bufferStatus = SR_OK;
    while ((readerStatus = SR_BamInStreamLoadPair(&pAlgnOne, &pAlgnTwo, pBamInStream)) == SR_OK)
    {
        SR_AlgnType algnType = SR_GetAlignmentType(&pAlgnOne, &pAlgnTwo, scTolerance, maxMismatchRate, minMQ);
        if ((algnType == SR_UNIQUE_ORPHAN || algnType == SR_UNIQUE_SOFT || algnType == SR_UNIQUE_MULTIPLE) && pBamInStream->numThreads > 0)
        {
            SR_BamInStreamPush(pBamInStream, pAlgnOne, threadID);
            bufferStatus = SR_BamInStreamPush(pBamInStream, pAlgnTwo, threadID);

            SR_BamInStreamSetAlgnType(pBamInStream, threadID, algnType);
        }
        else
        {
            SR_BamInStreamRecycle(pBamInStream, pAlgnOne);
            SR_BamInStreamRecycle(pBamInStream, pAlgnTwo);
        }

        if (bufferStatus == SR_FULL)
            break;
    }

    return readerStatus;
}

unsigned int SR_BamInStreamShrinkPool(SR_BamInStream* pBamInStream, unsigned int newSize)
{
    unsigned int currSize = pBamInStream->pMemPool->numBuffs;
    if (currSize > newSize)
    {
        khash_t(buffAddress)* buffHash = kh_init(buffAddress);
        kh_resize(buffAddress, buffHash, currSize);

        int ret = 0;
        khiter_t khIter = 0;
        int64_t address = 0;

        for (SR_BamNode* pUsedNode = pBamInStream->pAlgnLists[PREV_BIN].first; pUsedNode != NULL; pUsedNode = pUsedNode->next)
        {
            address = (int64_t) pUsedNode->whereFrom;
            kh_put(buffAddress, buffHash, address, &ret);
        }

        for (SR_BamNode* pUsedNode = pBamInStream->pAlgnLists[CURR_BIN].first; pUsedNode != NULL; pUsedNode = pUsedNode->next)
        {
            address = (int64_t) pUsedNode->whereFrom;
            kh_put(buffAddress, buffHash, address, &ret);
        }

        unsigned int delNum = currSize - newSize;
        SR_BamBuff* pPrevBuff = NULL;
        SR_BamBuff* pCurrBuff = pBamInStream->pMemPool->pFirstBuff;

        while (pCurrBuff != NULL && delNum != 0)
        {
            int64_t address = (int64_t) pCurrBuff;
            khIter = kh_get(buffAddress, buffHash, address);

            if (khIter == kh_end(buffHash))
            {
                SR_BamBuff* pDelBuff = pCurrBuff;
                pCurrBuff = pCurrBuff->nextBuff;

                SR_BamBuffClear(pDelBuff, pBamInStream->pMemPool);
                SR_BamBuffFree(pDelBuff, pBamInStream->pMemPool->buffCapacity);
                --delNum;
                --(pBamInStream->pMemPool->numBuffs);

                if (pPrevBuff != NULL)
                    pPrevBuff->nextBuff = pCurrBuff;
                else
                    pBamInStream->pMemPool->pFirstBuff = pCurrBuff;
            }
            else
            {
                pPrevBuff = pCurrBuff;
                pCurrBuff = pCurrBuff->nextBuff;
            }
        }

        kh_destroy(buffAddress, buffHash);
    }

    return pBamInStream->pMemPool->numBuffs;
}

