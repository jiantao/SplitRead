/*
 * =====================================================================================
 *
 *       Filename:  SR_FragLenDstrb.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/15/2011 03:43:08 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "khash.h"
#include "SR_Error.h"
#include "SR_FragLenDstrb.h"

// default capacity of fragment length distribution
#define DEFAULT_FRAG_DSTRB_CAP 10

// read group name hash
KHASH_MAP_INIT_STR(readGrpName, uint32_t);


static void SR_FragLenHistInit(SR_FragLenHist* pHist, uint32_t fragLen)
{
    pHist->min = fragLen;
    pHist->max = fragLen;

    pHist->lowerBound = fragLen / 2;
    pHist->size = fragLen * 2 - pHist->lowerBound;

    pHist->bin = (uint64_t*) calloc(pHist->size, sizeof(uint64_t));
    if (pHist->bin == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the fragment length histogram object.\n");
}

static void SR_FragLenHistResize(SR_FragLenHist* pHist, uint32_t newBound)
{
    if (newBound >= pHist->lowerBound + pHist->size)
    {
        uint32_t newSize = newBound - pHist->lowerBound;
        pHist->bin = (uint64_t*) realloc(pHist->bin, newSize * sizeof(uint64_t));
        if (pHist->bin == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of bins in the histogram object.\n");

        memset(pHist->bin + pHist->size, 0, newSize - pHist->size);
        pHist->size = newSize;
    }
    else if (newBound < pHist->lowerBound)
    {
        size_t newSize = pHist->lowerBound + pHist->size - newBound;
        uint64_t* newBin = calloc(newSize, sizeof(uint64_t));
        if (newBin == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of bins in the histogram object.\n");

        memcpy(newBin + (newSize - pHist->size), pHist->bin, pHist->size * sizeof(uint64_t));
        free(pHist->bin);
        pHist->bin = newBin;

        pHist->lowerBound = newBound;
        pHist->size = newSize;
    }
}

SR_FragLenDstrb* SR_FragLenDstrbAlloc(unsigned short minMQ, uint32_t capacity)
{
    SR_FragLenDstrb* pNewDstrb = (SR_FragLenDstrb*) calloc(1, sizeof(*pNewDstrb));
    if (pNewDstrb == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the fragment length distribution object.\n");

    if (capacity == 0)
        capacity = DEFAULT_FRAG_DSTRB_CAP;

    pNewDstrb->pReadGrpNames = (char**) calloc(capacity, sizeof(char*));
    if (pNewDstrb->pReadGrpNames == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of read group names in the fragment length distribution object.\n");


    pNewDstrb->pReadGrpHash = NULL;
    pNewDstrb->pHists = NULL;

    pNewDstrb->size = 0;
    pNewDstrb->capacity = capacity;

    pNewDstrb->minMQ = minMQ;
    pNewDstrb->hasRG = TRUE;

    return pNewDstrb;
}

void SR_FragLenDstrbFree(SR_FragLenDstrb* pDstrb)
{
    if (pDstrb != NULL)
    {
        free(pDstrb->pHists);
        if (pDstrb->pHists != NULL)
        {
            for (unsigned int i = 0; i != pDstrb->size; ++i)
                free(pDstrb->pHists[i].bin);

            free(pDstrb->pHists);
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

SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader)
{
    const char* readGrpPos = pBamHeader->pOrigHeader->text;
    while ((readGrpPos = strstr(readGrpPos, "@RG")) != NULL)
    {
        // expand the array if necessary
        if (pDstrb->size == pDstrb->capacity)
        {
            pDstrb->capacity *= 2;
            pDstrb->pReadGrpNames = (char**) realloc(pDstrb->pReadGrpNames, pDstrb->capacity * sizeof(char*));
            if (pDstrb->pReadGrpNames == NULL)
                SR_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");
        }

        // get the name of the current read group
        const char* readGrpNamePos = strstr(readGrpPos, "ID:");
        if (readGrpNamePos != NULL)
        {
            readGrpNamePos += 3;
        }
        else
        {
            SR_ErrMsg("ERROR: the \"ID\" field is required under the read group tag.\n");
            return SR_ERR;
        }

        const char* readGrpNameEnd = strpbrk(readGrpPos, " \t\n\0");
        size_t readGrpNameLen = readGrpNameEnd - readGrpNamePos;
        if (readGrpNameLen == 0)
        {
            SR_ErrMsg("ERROR: the length of read group ID is zero.\n");
            return SR_ERR;
        }

        pDstrb->pReadGrpNames[pDstrb->size] = (char*) calloc(readGrpNameLen + 1, sizeof(char));
        if (pDstrb->pReadGrpNames[pDstrb->size] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");

        strncpy(pDstrb->pReadGrpNames[pDstrb->size], readGrpPos, readGrpNameLen);

        ++(pDstrb->size);
        ++readGrpPos;
    }

    if (pDstrb->size != 0) // read group names are found. we insert them into a hash table
    {
        pDstrb->hasRG = TRUE;
        pDstrb->pReadGrpHash = (khash_t(readGrpName)*) kh_init(readGrpName);
        kh_resize(readGrpName, pDstrb->pReadGrpHash, pDstrb->capacity);

        int khRet = 0;
        khiter_t khIter;
        for (unsigned int i = 0; i != pDstrb->size; ++i)
        {
            khIter = kh_put(readGrpName, pDstrb->pReadGrpHash, pDstrb->pReadGrpNames[i], &khRet);
            if (khRet == 0)
            {
                SR_ErrMsg("ERROR: Found a non-unique read group ID.\n");
                return SR_ERR;
            }

            kh_value((khash_t(readGrpName)*) pDstrb->pReadGrpHash, khIter) = i;
        }
    }
    else // did not find any read group names. treat all the alignments as one read group
    {
        SR_ErrMsg("WARNING: No read group is found. The alignments in this bam file will be treated as from one read group.\n");

        pDstrb->size = 1;
        pDstrb->hasRG = FALSE;
    }

    pDstrb->pHists = (SR_FragLenHist*) calloc(pDstrb->size,  sizeof(SR_FragLenHist));
    if (pDstrb->pHists == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storageof the histogram in the fragment length distribution object.\n");

    return SR_OK;
}

SR_Status SR_FragLenDstrbUpdate(SR_FragLenDstrb* pDstrb, const SR_BamPairStats* pPairStats)
{
    unsigned int dstrbIndex = 0;

    if (pDstrb->hasRG)
    {
        khash_t(readGrpName)* pRgHash = pDstrb->pReadGrpHash;
        khiter_t khIter = kh_get(readGrpName, pRgHash, pPairStats->RG);
        if (khIter != kh_end(pRgHash))
        {
            dstrbIndex = kh_value(pRgHash, khIter);
        }
        else
        {
            SR_ErrMsg("ERROR: Found a read group name that is not record in the header\n");
            return SR_ERR;
        }

        pDstrb->pReadGrpHash = pRgHash;
    }

    SR_FragLenHist* pCurrHist = pDstrb->pHists + dstrbIndex;
    if (pCurrHist->bin == NULL)
    {
        SR_FragLenHistInit(pCurrHist, pPairStats->fragLen);
    }

    if (pPairStats->fragLen >= pCurrHist->lowerBound + pCurrHist->size)
    {
        SR_FragLenHistResize(pCurrHist, pPairStats->fragLen * 2);
    }
    else if (pPairStats->fragLen < pCurrHist->lowerBound)
    {
        SR_FragLenHistResize(pCurrHist, pPairStats->fragLen / 2);
    }

    ++(pCurrHist->bin[pPairStats->fragLen - pCurrHist->lowerBound]);
    ++(pCurrHist->total);
    ++(pDstrb->pairModeCount[SR_PairModeMap[pPairStats->pairMode]]);

    if (pPairStats->fragLen > pCurrHist->max)
    {
        pCurrHist->max = pPairStats->fragLen;
    }
    else if (pPairStats->fragLen < pCurrHist->min)
    {
        pCurrHist->min = pPairStats->fragLen;
    }

    if (pCurrHist->bin[pPairStats->fragLen - pCurrHist->lowerBound] > pCurrHist->bin[pCurrHist->mode])
    {
        pCurrHist->mode = pPairStats->fragLen - pCurrHist->lowerBound;
    }

    return SR_OK;
}
