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

#include <assert.h>

#include "khash.h"
#include "SR_Error.h"
#include "SR_Utilities.h"
#include "SR_QueryRegion.h"
#include "SR_BamInStream.h"


//===============================
// Type and constant definition
//===============================

// value indicate that we did not load any read
// into the bam array
#define NO_QUERY_YET (-2)

// default capacity of a bam array
#define DEFAULT_BAM_ARRAY_CAP 200

// default soft clipping tolerance
#define DEFAULT_SC_TOLERANCE 0.2

// default capacity of fragment length distribution
#define DEFAULT_FRAG_DSTRB_CAP 10

// alignment status
enum AlignmentStatus
{
    NEITHER_GOOD = -1,    // neither a good anchor nor a good orphan candidate

    GOOD_ANCHOR  = 0,     // a good anchor candidate

    GOOD_ORPHAN  = 1      // a good orphan candidate
};

// a mask used to filter out those unwanted reads for split alignments
// it includes proper paired reads, secondar reads, qc-failed reads and duplicated reads
#define SR_BAM_FMASK (BAM_FPROPER_PAIR | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)

#define PREV_BIN 0

#define CURR_BIN 1

typedef struct SR_BamNode SR_BamNode;

typedef struct SR_BamList SR_BamList;

typedef struct SR_BamBuff SR_BamBuff;

typedef struct SR_BamMemPool SR_BamMemPool;

struct SR_BamNode
{
    bam1_t alignment;

    SR_BamNode* prev;

    SR_BamNode* next;

    SR_BamBuff* whereFrom;
};

struct SR_BamList
{
    unsigned int numNode;

    SR_BamNode* first;

    SR_BamNode* last;
};

struct SR_BamBuff
{
    unsigned int numUsed;

    SR_BamBuff* nextBuff;

    SR_BamNode* pNodeArray;
};

struct SR_BamMemPool
{
    unsigned short numBuffs;

    unsigned int buffCapacity;

    SR_BamBuff* pFirstBuff;

    SR_BamNode* pFirstAvlNode;
};

// array of bam alignments
typedef struct SR_BamArray
{
    bam1_t* data;            // a pointer to the read-in bam alignments

    unsigned int size;       // how many alignments are in the array

    unsigned int capacity;   // maximum number of alignments the array can hold

}SR_BamArray;

// initialize a string-to-bam hash table used to retrieve a pair of read
KHASH_MAP_INIT_STR(queryName, SR_BamNode*);

KHASH_MAP_INIT_STR(readGrpName, uint32_t);

// private data structure that holds all bam-input-related information
struct SR_BamInStreamPrvt
{
    bamFile fpBamInput;                        // file pointer to a input bam file

    bam_index_t* pBamIndex;                    // file pointer to a input bam index file

    SR_BamMemPool* pMemPool;                   // memory pool used to allocate and recycle the bam alignments

    khash_t(queryName)* pNameHashes[2];     

    SR_BamList* pRetLists;

    SR_BamNode* pNewNode;

    SR_BamList pAlgnLists[2];                  // lists used to store those incoming alignments. each thread has their own lists.

    unsigned int numThreads;                   // number of threads will be used

    unsigned int reportSize;                   // number of pairs should be loaded before report

    int32_t currRefID;                         // the reference ID of the current read-in alignment

    int32_t currBinPos;                        // the start position of current bin (0-based)

    uint32_t binLen;                           // the length of bin

    double scTolerance;                        // soft clipping tolerance rate. 
};


//===================
// Static functions
//===================

/* 

static SR_BamArray* SR_BamArrayAlloc(unsigned int capacity)
{
    SR_BamArray* pBamArray = (SR_BamArray*) malloc(sizeof(SR_BamArray));
    if (pBamArray == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a bam array object");

    pBamArray->data = (bam1_t*) calloc(capacity, sizeof(bam1_t));
    if (pBamArray->data == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storeage of bam alignments in a bam array object");

    pBamArray->size = 0;
    pBamArray->capacity = capacity;

    return pBamArray;
}

static void SR_BamArrayFree(SR_BamArray* pBamArray)
{
    if (pBamArray != NULL)
    {
        if (pBamArray->data != NULL)
        {
            for (unsigned int i = 0; i != pBamArray->capacity; ++i)
                free(pBamArray->data[i].data);

            free(pBamArray->data);
        }

        free(pBamArray);
    }
}

static void SR_BamArrayStartOver(SR_BamArray* pBamArray)
{
    if (SR_ARRAY_GET_SIZE(pBamArray) <= 1)
        return;

    bam_copy1(SR_ARRAY_GET_PT(pBamArray, 0), SR_ARRAY_GET_LAST_PT(pBamArray));
    pBamArray->size = 1;
}

*/

static SR_BamBuff* SR_BamBuffAlloc(unsigned int buffCapacity)
{
    SR_BamBuff* pNewBuff = (SR_BamBuff*) malloc(sizeof(SR_BamBuff));
    if (pNewBuff == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the bam memory node object.\n");

    pNewBuff->numUsed = 0;
    pNewBuff->nextBuff = NULL;

    pNewBuff->pNodeArray = (SR_BamNode*) calloc(buffCapacity, sizeof(SR_BamNode));
    if (pNewBuff->pNodeArray == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of alignments in the bam memory node object.\n");

    for (unsigned int i = 0; i != buffCapacity - 1; ++i)
    {
        pNewBuff->pNodeArray[i].whereFrom = pNewBuff;
        pNewBuff->pNodeArray[i].next = &(pNewBuff->pNodeArray[i + 1]);
    }

    pNewBuff->pNodeArray[buffCapacity - 1].whereFrom = pNewBuff;

    return pNewBuff;
}

static void SR_BamBuffFree(SR_BamBuff* pBuff)
{
    if (pBuff != NULL)
    {
        free(pBuff->pNodeArray);
        free(pBuff);
    }
}

static void SR_BamMemPoolExpand(SR_BamMemPool* pMemPool)
{
    SR_BamBuff* pNewBuff = SR_BamBuffAlloc(pMemPool->buffCapacity);

    pNewBuff->nextBuff = pMemPool->pFirstBuff;
    pMemPool->pFirstBuff = pNewBuff;

    pNewBuff->pNodeArray[pMemPool->buffCapacity - 1].next = pMemPool->pFirstAvlNode;
    pMemPool->pFirstAvlNode = &(pNewBuff->pNodeArray[0]);

    ++(pMemPool->numBuffs);
}

static SR_BamMemPool* SR_BamMemPoolAlloc(unsigned int buffCapacity)
{
    SR_BamMemPool* pNewPool = (SR_BamMemPool*) malloc(sizeof(SR_BamMemPool));
    if (pNewPool == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the bam memory pool object.\n");

    pNewPool->numBuffs = 0;
    pNewPool->pFirstAvlNode = NULL;
    pNewPool->pFirstBuff = NULL;
    pNewPool->buffCapacity = buffCapacity;

    SR_BamMemPoolExpand(pNewPool);

    return pNewPool;
}

static void SR_BamMemPoolFree(SR_BamMemPool* pMemPool)
{
    if (pMemPool != NULL)
    {
        SR_BamBuff* pDelBuff = pMemPool->pFirstBuff;
        SR_BamBuff* pNextBuff = pMemPool->pFirstBuff;
        while (pDelBuff != NULL)
        {
            pNextBuff = pDelBuff->nextBuff;
            SR_BamBuffFree(pDelBuff);
            pDelBuff = pNextBuff;
        }

        free(pMemPool);
    }
}

static void SR_BamListPushHead(SR_BamList* pList, SR_BamNode* pNewFirstNode)
{
    pNewFirstNode->next = pList->first;
    pNewFirstNode->prev = NULL;

    if (pList->first != NULL)
        pList->first->prev = pNewFirstNode;
    else
        pList->last = pNewFirstNode;

    pList->first = pNewFirstNode;
    ++(pList->numNode);
}

static void SR_BamListRemove(SR_BamList* pList, SR_BamNode* pNode)
{
    if (pNode->next != NULL)
        pNode->next->prev = pNode->prev;
    else
        pList->last = pNode->prev;

    if (pNode->prev != NULL)
        pNode->prev->next = pNode->next;
    else
        pList->first = pNode->next;

    pNode->next = NULL;
    pNode->prev = NULL;

    --(pList->numNode);
}

static void SR_BamListReset(SR_BamList* pList, SR_BamMemPool* pMemPool)
{
    if (pList->numNode == 0)
        return;

    pList->last->next = pMemPool->pFirstAvlNode;
    pMemPool->pFirstAvlNode = pList->first;

    pList->first = NULL;
    pList->last = NULL;
    pList->numNode = 0;
}

static SR_BamNode* SR_BamNodeAlloc(SR_BamMemPool* pMemPool)
{
    if (pMemPool->pFirstAvlNode == NULL)
    {
        SR_BamMemPoolExpand(pMemPool);
    }

    SR_BamNode* pNewNode = pMemPool->pFirstAvlNode;
    pMemPool->pFirstAvlNode = pMemPool->pFirstAvlNode->next;

    pNewNode->prev = NULL;
    pNewNode->next = NULL;

    return pNewNode;
}

static void SR_BamNodeFree(SR_BamNode* pNode, SR_BamMemPool* pMemPool)
{
    pNode->prev = NULL;
    pNode->next = pMemPool->pFirstAvlNode;
    pMemPool->pFirstAvlNode = pNode;
}

static inline int SR_BamInStreamLoadNext(SR_BamInStream* pBamInStream)
{
    // for the bam alignment array, if we need to expand its space
    // we have to initialize those newly created bam alignment 
    // and update the query name hash since the address of those
    // bam alignments are changed after expanding
    pBamInStream->pNewNode = SR_BamNodeAlloc(pBamInStream->pMemPool);
    int ret = bam_read1(pBamInStream->fpBamInput, &(pBamInStream->pNewNode->alignment));

    return ret;
}

static int SR_CheckSC(bam1_t* pAlignment, double scTolerance)
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
        if (pAlignment->core.qual != 0)
            return GOOD_ANCHOR;
        else
            return NEITHER_GOOD;
    }
    else if (isHeadSC && isTailSC)
        return NEITHER_GOOD;
    else
        return GOOD_ORPHAN;

}

static SR_Bool SR_IsQualifiedPair(SR_BamNode** ppAnchor, SR_BamNode** ppOrphan, double scTolerance)
{
    int anchorStatus = SR_CheckSC(&((*ppAnchor)->alignment), scTolerance);
    int orphanStatus = SR_CheckSC(&((*ppOrphan)->alignment), scTolerance);

    if ((anchorStatus == NEITHER_GOOD || orphanStatus == NEITHER_GOOD)
        || (anchorStatus == GOOD_ANCHOR && orphanStatus == GOOD_ANCHOR)
        || (anchorStatus == GOOD_ORPHAN && orphanStatus == GOOD_ORPHAN))
    {
        return FALSE;
    }
    else if (anchorStatus == GOOD_ORPHAN)
        SR_SWAP(*ppAnchor, *ppOrphan, SR_BamNode*);

    return TRUE;
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


//===============================
// Constructors and Destructors
//===============================

SR_BamInStream* SR_BamInStreamAlloc(const char* bamFilename, uint32_t binLen, unsigned int numThreads, unsigned int buffCapacity, unsigned int reportSize, double scTolerance)
{
    SR_BamInStream* pBamInStream = (SR_BamInStream*) calloc(1, sizeof(struct SR_BamInStreamPrvt));
    if (pBamInStream == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a bam input stream object.");

    pBamInStream->fpBamInput = bam_open(bamFilename, "r");
    if (pBamInStream->fpBamInput == NULL)
        SR_ErrQuit("ERROR: Cannot open bam file %s for reading.\n", bamFilename);

    pBamInStream->pBamIndex = bam_index_load(bamFilename);
    if (pBamInStream->pBamIndex == NULL)
        SR_ErrMsg("WARNING: Cannot open bam index file for reading. No jump allowed.\n");

    pBamInStream->numThreads = numThreads;
    pBamInStream->reportSize = reportSize;
    pBamInStream->currRefID = NO_QUERY_YET;
    pBamInStream->currBinPos = NO_QUERY_YET;
    pBamInStream->binLen = binLen;
    pBamInStream->scTolerance = scTolerance;
    pBamInStream->pNewNode = NULL;

    pBamInStream->pRetLists = (SR_BamList*) calloc(numThreads, sizeof(SR_BamList));
    if (pBamInStream->pRetLists == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of retrun alignment lists in the bam input stream object.\n");

    pBamInStream->pNameHashes[PREV_BIN] = kh_init(queryName);
    kh_resize(queryName, pBamInStream->pNameHashes[PREV_BIN], reportSize);

    pBamInStream->pNameHashes[CURR_BIN] = kh_init(queryName);
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
        SR_BamMemPoolFree(pBamInStream->pMemPool);

        bam_close(pBamInStream->fpBamInput);
        bam_index_destroy(pBamInStream->pBamIndex);

        free(pBamInStream);
    }
}

SR_BamHeader* SR_BamHeaderAlloc(void)
{
    SR_BamHeader* pNewHeader = (SR_BamHeader*) calloc(1, sizeof(SR_BamHeader));
    if (pNewHeader == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a bam header object");

    return pNewHeader;
}

void SR_BamHeaderFree(SR_BamHeader* pBamHeader)
{
    if (pBamHeader != NULL)
    {
        free(pBamHeader->pMD5s);
        bam_header_destroy(pBamHeader->pOrigHeader);

        free(pBamHeader);
    }
}

SR_FragLenDstrb* SR_FragLenDstrbAlloc(uint8_t drctField, uint32_t capacity)
{
    SR_FragLenDstrb* pNewDstrb = (SR_FragLenDstrb*) malloc(sizeof(*pNewDstrb));
    if (pNewDstrb == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the fragment length distribution object.\n");

    if (capacity == 0)
        capacity = DEFAULT_FRAG_DSTRB_CAP;

    pNewDstrb->pReadGrpNames = (char**) calloc(capacity, sizeof(char*));
    if (pNewDstrb->pReadGrpNames == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of read group names in the fragment length distribution object.\n");

    pNewDstrb->pReadGrpHash = kh_init(readGrpName);
    pNewDstrb->pFragLenHists = NULL;

    pNewDstrb->drctField = drctField;
    pNewDstrb->size = 0;
    pNewDstrb->capacity = capacity;

    return pNewDstrb;
}

void SR_FragLenDstrbFree(SR_FragLenDstrb* pDstrb)
{
    if (pDstrb != NULL)
    {
        if (pDstrb->pFragLenHists != NULL)
        {
            for (unsigned int i = 0; i != pDstrb->size; ++i)
                gsl_histogram_free(pDstrb->pFragLenHists[i]);

            free(pDstrb->pFragLenHists);
        }

        if (pDstrb->pReadGrpNames != NULL)
        {
            for (unsigned int i = 0; i != pDstrb->size; ++i)
                free(pDstrb->pReadGrpNames[i]);

            free(pDstrb->pReadGrpNames);
        }

        kh_destroy(readGrpName, pDstrb->pReadGrpHash);

        free(pDstrb);
    }
}

//======================
// Interface functions
//======================

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
                kh_value(pBamInStream->pNameHashes[CURR_BIN], khIter) = pBamInStream->pNewNode;
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

SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader)
{
    for (const char* readGrpPos = pBamHeader->pOrigHeader->text; ; ++readGrpPos)
    {
        if ((readGrpPos = strstr(readGrpPos, "@RG")) != NULL 
            && (readGrpPos = strstr(readGrpPos, "ID:")) != NULL)
        {
            const char* readGrpEnd = strpbrk(readGrpPos, " \t\n\0");
            if (readGrpEnd == NULL)
                return SR_ERR;

            if (pDstrb->size == pDstrb->capacity)
            {
                pDstrb->capacity *= 2;
                pDstrb->pReadGrpNames = (char**) realloc(pDstrb->pReadGrpNames, pDstrb->capacity * sizeof(char*));
                if (pDstrb->pReadGrpNames == NULL)
                    SR_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");
            }

            size_t readGrpLen = readGrpEnd - readGrpPos + 1;
            pDstrb->pReadGrpNames[pDstrb->size] = (char*) calloc(readGrpLen + 1, sizeof(char));
            if (pDstrb->pReadGrpNames[pDstrb->size] == NULL)
                SR_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");

            strncpy(pDstrb->pReadGrpNames[pDstrb->size], readGrpPos, readGrpLen);
            ++(pDstrb->size);
        }
        else
            break;
    }

    if (pDstrb->size == 0)
    {
        SR_FragLenDstrbFree(pDstrb);
        return SR_NOT_FOUND;
    }
    else
    {
        int khRet = 0;
        khiter_t khIter;
        for (unsigned int i = 0; i != pDstrb->size; ++i)
        {
            khIter = kh_put(readGrpName, pDstrb->pReadGrpHash, pDstrb->pReadGrpNames[i], &khRet);
            if (khRet == 0)
            {
                SR_ErrMsg("ERROR: Found a non-unique read group name.\n");
                return SR_ERR;
            }

            kh_value((khash_t(readGrpName)*) pDstrb->pReadGrpHash, khIter) = i;
        }
    }

    return SR_OK;
}

void SR_FragLenDstrbSetHist(SR_FragLenDstrb* pDstrb, size_t numBins)
{
    pDstrb->pFragLenHists = (gsl_histogram**) malloc(pDstrb->size * sizeof(gsl_histogram*));
    if (pDstrb->pFragLenHists == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storageof the histogram in the fragment length distribution object.\n");

    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        pDstrb->pFragLenHists[i] = gsl_histogram_alloc(numBins);
        if (pDstrb->pFragLenHists[i] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storageof the histogram in the fragment length distribution object.\n");
    }
}

// read an alignment from a bam file
SR_Status SR_BamInStreamRead(bam1_t* pAlignment, SR_BamInStream* pBamInStream)
{
    int ret = bam_read1(pBamInStream->fpBamInput, pAlignment);

    if (ret > 0)
        return SR_OK;
    else if (ret == -1)
        return SR_EOF;
    else
        return SR_ERR;
}

// get the current reference ID
int32_t SR_BamInStreamGetRefID(const SR_BamInStream* pBamInStream)
{
    return (pBamInStream->currRefID);
}

// load a unique-orphan pair from a bam file
SR_Status SR_BamInStreamGetPairs(SR_BamInStream* pBamInStream, unsigned int threadID, SR_FragLenDstrb* pDstrb)
{
    SR_BamList* pLoadingList = pBamInStream->pRetLists + threadID;

    int ret = 1;
    while(ret > 0 && (ret = SR_BamInStreamLoadNext(pBamInStream)) > 0)
    {
        // exclude those reads who are non-paired-end, qc-fail, duplicate-marked, proper-paired, 
        // both aligned, secondary-alignment and no-name-specified.
        if ((pBamInStream->pNewNode->alignment.core.flag & SR_BAM_FMASK) != 0
            || strcmp(bam1_qname(&(pBamInStream->pNewNode->alignment)), "*") == 0)
        {
            SR_BamNodeFree(pBamInStream->pNewNode, pBamInStream->pMemPool);
            pBamInStream->pNewNode = NULL;
            continue;
        }

        // update the current ref ID or position if the incoming alignment has a 
        // different value. The name hash and the bam array will be reset
        if (pBamInStream->pNewNode->alignment.core.tid != pBamInStream->currRefID
            || pBamInStream->pNewNode->alignment.core.pos >= pBamInStream->currBinPos + 2 * pBamInStream->binLen)
        {
            pBamInStream->currRefID  = pBamInStream->pNewNode->alignment.core.tid;
            pBamInStream->currBinPos = pBamInStream->pNewNode->alignment.core.pos;

            kh_clear(queryName, pBamInStream->pNameHashes[PREV_BIN]);
            kh_clear(queryName, pBamInStream->pNameHashes[CURR_BIN]);

            SR_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);
            SR_BamListReset(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pMemPool);

            if (pBamInStream->pNewNode->alignment.core.tid != pBamInStream->currRefID)
            {
                ret = SR_OUT_OF_RANGE;
            }
        }
        else if (pBamInStream->pNewNode->alignment.core.pos >= pBamInStream->currBinPos + pBamInStream->binLen)
        {
            pBamInStream->currBinPos += pBamInStream->binLen;

            kh_clear(queryName, pBamInStream->pNameHashes[PREV_BIN]);
            SR_SWAP(pBamInStream->pNameHashes[PREV_BIN], pBamInStream->pNameHashes[CURR_BIN], khash_t(queryName)*);

            SR_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);

            SR_SWAP(pBamInStream->pAlgnLists[PREV_BIN], pBamInStream->pAlgnLists[PREV_BIN], SR_BamList);
        }

        SR_BamListPushHead(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pNewNode);

        SR_BamNode* pAnchor = NULL;
        SR_BamNode* pOrphan = NULL;

        khiter_t khIter;
        khIter = kh_get(queryName, pBamInStream->pNameHashes[PREV_BIN], bam1_qname(&(pBamInStream->pNewNode->alignment)));

        if (khIter != kh_end(pBamInStream->pNameHashes[PREV_BIN]))
        {
            pAnchor = kh_value(pBamInStream->pNameHashes[PREV_BIN], khIter);
            pOrphan = pBamInStream->pNewNode;

            kh_del(queryName, pBamInStream->pNameHashes[PREV_BIN], khIter);

            SR_BamListRemove(&(pBamInStream->pAlgnLists[PREV_BIN]), pAnchor);
            SR_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), pOrphan);

            if (SR_IsQualifiedPair(&pAnchor, &pOrphan, pBamInStream->scTolerance))
            {
                SR_BamListPushHead(pLoadingList, pOrphan);
                SR_BamListPushHead(pLoadingList, pAnchor);

                if (pLoadingList->numNode == pBamInStream->reportSize * 2)
                    ret = SR_OK;
            }
            else
            {
                SR_BamNodeFree(pAnchor, pBamInStream->pMemPool);
                SR_BamNodeFree(pOrphan, pBamInStream->pMemPool);
            }
        }
        else
        {
            int khRet = 0;
            khIter = kh_put(queryName, pBamInStream->pNameHashes[CURR_BIN], bam1_qname(&(pBamInStream->pNewNode->alignment)), &khRet);

            if (khRet == 0) // we found a pair of alignments 
            {
                pAnchor = kh_value(pBamInStream->pNameHashes[CURR_BIN], khIter);
                pOrphan = pBamInStream->pNewNode;

                kh_del(queryName, pBamInStream->pNameHashes[CURR_BIN], khIter);

                SR_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), pAnchor);
                SR_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), pOrphan);

                if (SR_IsQualifiedPair(&pAnchor, &pOrphan, pBamInStream->scTolerance))
                {
                    SR_BamListPushHead(pLoadingList, pOrphan);
                    SR_BamListPushHead(pLoadingList, pAnchor);

                    if (pLoadingList->numNode == pBamInStream->reportSize * 2)
                        ret = SR_OK;
                }
                else
                {
                    SR_BamNodeFree(pAnchor, pBamInStream->pMemPool);
                    SR_BamNodeFree(pOrphan, pBamInStream->pMemPool);
                }
            }
            else // not finding corresponding mate, save the current value and move on
            {
                kh_value(pBamInStream->pNameHashes[CURR_BIN], khIter) = pBamInStream->pNewNode;
            }
        }
    }

    if (ret < 0)
    {
        if ( ret != SR_OUT_OF_RANGE && ret != SR_EOF)
            return SR_ERR;
    }

    return ret;
}
