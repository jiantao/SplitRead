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

#include <stdio.h>
#include <math.h>

#include "khash.h"
#include "SR_Error.h"
#include "SR_FragLenDstrb.h"

// default capacity of fragment length distribution
#define DEFAULT_FRAG_DSTRB_CAP 10

typedef struct SR_FragLenBin
{
    uint32_t fragLen;

    uint64_t freq;

}SR_FragLenBin;


// read group name hash
KHASH_MAP_INIT_STR(readGrpName, uint32_t);

// fragment length hash
KHASH_MAP_INIT_INT(fragLen, uint64_t);

static inline void SR_FragLenHistClearRaw(SR_FragLenHist* pHist)
{
    for (unsigned int i = 0; i != NUM_PAIR_MODE; ++i)
    {
        kh_destroy(fragLen, pHist->rawHist[i]);
        pHist->rawHist[i] = NULL;
    }
}

static inline void SR_FragLenHistClear(SR_FragLenHist* pHist)
{
    SR_FragLenHistClearRaw(pHist);

    for (unsigned int i = 0; i != NUM_TOP_PAIR_MODE; ++i)
    {
        free(pHist->fragLen);
        free(pHist->cdf);

        pHist->fragLen = NULL;
        pHist->cdf = NULL;
    }
}

// we sort the pair mode count array decreasingly 
static int ComparePairModeBin(const void* a, const void* b)
{
    const SR_PairModeBin* pModeBinOne = a;
    const SR_PairModeBin* pModeBinTwo = b;

   if (pModeBinOne->freq > pModeBinTwo->freq)
       return -1;
   else if (pModeBinOne->freq < pModeBinTwo->freq)
       return 1;
   else
       return 0;
}

static int CompareFragLenBin(const void* a, const void* b)
{
    const SR_FragLenBin* pBinOne = a;
    const SR_FragLenBin* pBinTwo = b;

   if (pBinOne->fragLen > pBinTwo->fragLen)
       return 1;
   else if (pBinOne->fragLen < pBinTwo->fragLen)
       return -1;
   else
       return 0;
}

static inline SR_Bool SR_IsValidPairModeSet(const SR_PairModeBin* pBestBin, const SR_PairModeBin* pSecBestBin)
{
    unsigned int pairModeSet = (pBestBin->pairMode << 3 | pSecBestBin->pairMode);
    if (SR_PairModeSetMap[pairModeSet] == 1)
        return TRUE;
    else
        return FALSE;
}

static void SR_FragLenHistMergeRaw(SR_FragLenHist* pHist)
{
    khash_t(fragLen)* pBestRaw = pHist->rawHist[pHist->modeCount[0].pairMode];
    const khash_t(fragLen)* pSecBestRaw = pHist->rawHist[pHist->modeCount[1].pairMode];

    int ret = 0;
    khiter_t bestIter = 0;
    for (khiter_t secBestIter = kh_begin(pSecBestRaw); secBestIter != kh_end(pSecBestRaw); ++secBestIter)
    {
        if (kh_exist(pSecBestRaw, secBestIter))
        {
            bestIter = kh_put(fragLen, pBestRaw, kh_key(pSecBestRaw, secBestIter), &ret);
            if (ret == 0)
                kh_value(pBestRaw, bestIter) += kh_value(pSecBestRaw, secBestIter);
            else
                kh_value(pBestRaw, bestIter) = kh_value(pSecBestRaw, secBestIter);
        }
    }

    pHist->rawHist[pHist->modeCount[0].pairMode] = pBestRaw;
}

static void SR_FragLenHistToMature(SR_FragLenHist* pHist)
{
    const khash_t(fragLen)* pBestRaw = pHist->rawHist[pHist->modeCount[0].pairMode];
    SR_FragLenBin* matureHist = (SR_FragLenBin*) malloc(kh_size(pBestRaw) * sizeof(SR_FragLenBin));
    if (matureHist == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the histogram in the fragment length distribution object.\n");

    pHist->size = kh_size(pBestRaw);

    unsigned int i = 0;
    for (khiter_t khIter = kh_begin(pBestRaw); khIter != kh_end(pBestRaw); ++khIter)
    {
        if (kh_exist(pBestRaw, khIter))
        {
            matureHist[i].fragLen = kh_key(pBestRaw, khIter);
            matureHist[i].freq = kh_value(pBestRaw, khIter);
            ++i;
        }
    }

    qsort(matureHist, pHist->size, sizeof(SR_FragLenBin), CompareFragLenBin);

    pHist->fragLen = (uint32_t*) malloc(pHist->size * sizeof(uint32_t));
    if (pHist->fragLen == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storeage of the fragment length array in the fragment length histogram.\n");

    pHist->cdf = (double*) malloc(pHist->size * sizeof(double));
    if (pHist->cdf == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storeage of the cdf in the fragment length histogram.\n");

    double cumFreq = 0.0;
    double totalFragLen = 0.0;
    uint64_t totalFreq = pHist->modeCount[0].freq + pHist->modeCount[1].freq;
    for (unsigned int j = 0; j != pHist->size; ++j)
    {
        totalFragLen += matureHist[j].fragLen * matureHist[j].freq;
        pHist->fragLen[j] = matureHist[j].fragLen;
        cumFreq += matureHist[j].freq;
        pHist->cdf[j] = cumFreq / totalFreq;
    }

    pHist->mean = totalFragLen / totalFreq;

    for (unsigned int j = 0; j != pHist->size; ++j)
    {
        if (pHist->cdf[j] >= 0.5)
        {
            pHist->median = pHist->fragLen[j];
            break;
        }
    }

    pHist->stdev = 0.0;
    for (unsigned int j = 0; j != pHist->size; ++j)
    {
        pHist->stdev += (double) matureHist[j].freq * pow(pHist->mean - matureHist[j].fragLen, 2);
    }

    if (totalFreq != 1)
        pHist->stdev = sqrt(pHist->stdev / (double) (totalFreq - 1));

    free(matureHist);
}


static SR_Status SR_FragLenHistSetMature(SR_FragLenHist* pHist)
{
    qsort(pHist->modeCount, NUM_PAIR_MODE, sizeof(SR_PairModeBin), ComparePairModeBin);

    SR_Bool isValid = SR_IsValidPairModeSet(&(pHist->modeCount[0]), &(pHist->modeCount[1]));

    if (isValid)
    {
        SR_FragLenHistMergeRaw(pHist);
        SR_FragLenHistToMature(pHist);
    }

    SR_FragLenHistClearRaw(pHist);

    if (isValid)
        return SR_OK;
    else
        return SR_ERR;
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
        if (pDstrb->pHists != NULL)
        {
            for (unsigned int i = 0; i != pDstrb->size; ++i) 
                SR_FragLenHistClear(&(pDstrb->pHists[i]));

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

        const char* readGrpNameEnd = strpbrk(readGrpNamePos, " \t\n\0");
        size_t readGrpNameLen = readGrpNameEnd - readGrpNamePos;
        if (readGrpNameLen == 0)
        {
            SR_ErrMsg("ERROR: the length of read group ID is zero.\n");
            return SR_ERR;
        }

        pDstrb->pReadGrpNames[pDstrb->size] = (char*) calloc(readGrpNameLen + 1, sizeof(char));
        if (pDstrb->pReadGrpNames[pDstrb->size] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");


        strncpy(pDstrb->pReadGrpNames[pDstrb->size], readGrpNamePos, readGrpNameLen);

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

    // initialize the pair mode in each histogram
    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        for (unsigned int j = 0; j != NUM_PAIR_MODE; ++j)
            pDstrb->pHists[i].modeCount[j].pairMode = j;
    }

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
    khash_t(fragLen)* pCurrHash = pCurrHist->rawHist[pPairStats->pairMode];

    if (pCurrHash == NULL)
    {
        pCurrHash = kh_init(fragLen);
        kh_resize(fragLen, pCurrHash, 20);
    }

    int ret = 0;
    khiter_t khIter = kh_put(fragLen, pCurrHash, pPairStats->fragLen, &ret);
    if (ret == 0)
    {
        kh_value(pCurrHash, khIter) += 1;
    }
    else
    {
        kh_value(pCurrHash, khIter) = 1;
    }

    ++(pCurrHist->modeCount[pPairStats->pairMode].freq);
    
    pCurrHist->rawHist[pPairStats->pairMode] = pCurrHash;

    return SR_OK;
}

void SR_FragLenDstrbFinalize(SR_FragLenDstrb* pDstrb)
{
    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        SR_Status pairSetStatus = SR_FragLenHistSetMature(&(pDstrb->pHists[i]));

        if (pairSetStatus == SR_ERR)
        {
            SR_ErrMsg("WARNING: Read group \"%s\" has incompatible pair mode set(%d and %d).\n"
                      "         Fragment length histogram will not be calculated for this group.\n", 
                       pDstrb->pReadGrpNames[i], pDstrb->pHists[i].modeCount[0].pairMode + 1, pDstrb->pHists[i].modeCount[1].pairMode + 1);
        }
    }
}

void SR_FragLenDstrbWrite(const SR_FragLenDstrb* pDstrb, FILE* dstrbOutput)
{
    size_t writeSize = 0;
    writeSize = fwrite(&(pDstrb->size), sizeof(uint32_t), 1, dstrbOutput);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");
    
    uint8_t rgFlag = pDstrb->hasRG;
    writeSize = fwrite(&(rgFlag), sizeof(uint8_t), 1, dstrbOutput);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");

    if (pDstrb->hasRG)
    {
        for (unsigned int i = 0; i != pDstrb->size; ++i)
        {
            uint16_t nameLen = strlen(pDstrb->pReadGrpNames[i]);
            writeSize = fwrite(&(nameLen), sizeof(uint16_t), 1, dstrbOutput);
            if (writeSize != 1)
                SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");
        }

        for (unsigned int i = 0; i != pDstrb->size; ++i)
        {
            uint32_t nameLen = strlen(pDstrb->pReadGrpNames[i]);
            writeSize = fwrite(pDstrb->pReadGrpNames[i], sizeof(char), nameLen, dstrbOutput);
            if (writeSize != nameLen)
                SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");
        }
    }

    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        uint8_t pairMode = pDstrb->pHists[i].modeCount[0].pairMode;
        writeSize = fwrite(&(pairMode), sizeof(uint8_t), 1, dstrbOutput);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");

        pairMode = pDstrb->pHists[i].modeCount[1].pairMode;
        writeSize = fwrite(&(pairMode), sizeof(uint8_t), 1, dstrbOutput);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");
    }

    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        uint64_t histFreq = pDstrb->pHists[i].modeCount[0].freq + pDstrb->pHists[i].modeCount[1].freq;
        writeSize = fwrite(&(histFreq), sizeof(uint64_t), 1, dstrbOutput);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");
    }

    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        writeSize = fwrite(&(pDstrb->pHists[i].mean), sizeof(double), 1, dstrbOutput);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");
    }

    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        writeSize = fwrite(&(pDstrb->pHists[i].median), sizeof(double), 1, dstrbOutput);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");
    }

    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        writeSize = fwrite(&(pDstrb->pHists[i].stdev), sizeof(double), 1, dstrbOutput);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");
    }

    uint32_t startIndex = 0;
    for (unsigned int i = 0; i != pDstrb->size + 1; ++i)
    {
        writeSize = fwrite(&(startIndex), sizeof(uint32_t), 1, dstrbOutput);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");

        startIndex += pDstrb->pHists[i].size;
    }

    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        writeSize = fwrite(pDstrb->pHists[i].fragLen, sizeof(uint32_t), pDstrb->pHists[i].size, dstrbOutput);
        if (writeSize != pDstrb->pHists[i].size)
            SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");
    }

    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        writeSize = fwrite(pDstrb->pHists[i].cdf, sizeof(double), pDstrb->pHists[i].size, dstrbOutput);
        if (writeSize != pDstrb->pHists[i].size)
            SR_ErrQuit("ERROR: Found an error when writing into the fragment length distribution file.\n");
    }

    fflush(dstrbOutput);
}
