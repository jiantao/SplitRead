/*
 * =====================================================================================
 *
 *       Filename:  SR_FragLenDstrb.c
 *
 *    jescription:  
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

/* 
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

static inline SR_Bool SR_IsValidPairModeSet(SR_PairMode modeOne, SR_PairMode modeTwo)
{
    unsigned int pairModeSet = (modeOne << 3 | modeTwo);
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
*/

static inline void SR_FragLenHistClearRaw(SR_FragLenHist* pHist)
{
    for (unsigned int i = 0; i != NUM_ALLOWED_HIST; ++i)
    {
        kh_destroy(fragLen, pHist->rawHist[i]);
        pHist->rawHist[i] = NULL;
    }
}

static inline void SR_FragLenHistClear(SR_FragLenHist* pHist)
{
    SR_FragLenHistClearRaw(pHist);

    for (unsigned int i = 0; i != NUM_ALLOWED_HIST; ++i)
    {
        free(pHist->fragLen[i]);
        free(pHist->cdf[i]);

        pHist->fragLen[i] = NULL;
        pHist->cdf[i] = NULL;
    }
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

static void SR_FragLenHistToMature(SR_FragLenHist* pHist)
{
    for (unsigned int k = 0; k != NUM_ALLOWED_HIST; ++k)
    {
        const khash_t(fragLen)* pRawHist = pHist->rawHist[k];
        if (pRawHist == NULL)
            continue;

        SR_FragLenBin* matureHist = (SR_FragLenBin*) malloc(kh_size(pRawHist) * sizeof(SR_FragLenBin));
        if (matureHist == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the histogram in the fragment length distribution object.\n");

        pHist->size[k] = kh_size(pRawHist);

        unsigned int i = 0;
        for (khiter_t khIter = kh_begin(pRawHist); khIter != kh_end(pRawHist); ++khIter)
        {
            if (kh_exist(pRawHist, khIter))
            {
                matureHist[i].fragLen = kh_key(pRawHist, khIter);
                matureHist[i].freq = kh_value(pRawHist, khIter);
                ++i;
            }
        }

        qsort(matureHist, pHist->size[k], sizeof(SR_FragLenBin), CompareFragLenBin);

        pHist->fragLen[k] = (uint32_t*) malloc(pHist->size[k] * sizeof(uint32_t));
        if (pHist->fragLen == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storeage of the fragment length array in the fragment length histogram.\n");

        pHist->cdf[k] = (double*) malloc(pHist->size[k] * sizeof(double));
        if (pHist->cdf == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storeage of the cdf in the fragment length histogram.\n");

        double cumFreq = 0.0;
        double totalFragLen = 0.0;
        uint64_t totalFreq = pHist->modeCount[k];
        for (unsigned int j = 0; j != pHist->size[k]; ++j)
        {
            totalFragLen += matureHist[j].fragLen * matureHist[j].freq;
            pHist->fragLen[k][j] = matureHist[j].fragLen;
            cumFreq += matureHist[j].freq;
            pHist->cdf[k][j] = cumFreq / totalFreq;
        }

        pHist->mean[k] = totalFragLen / totalFreq;

        for (unsigned int j = 0; j != pHist->size[k]; ++j)
        {
            if (pHist->cdf[k][j] >= 0.5)
            {
                pHist->median[k] = pHist->fragLen[k][j];
                break;
            }
        }

        pHist->stdev[k] = 0.0;
        for (unsigned int j = 0; j != pHist->size[k]; ++j)
        {
            pHist->stdev[k] += (double) matureHist[j].freq * pow(pHist->mean[k] - matureHist[j].fragLen, 2);
        }

        if (totalFreq != 1)
            pHist->stdev[k] = sqrt(pHist->stdev[k] / (double) (totalFreq - 1));

        free(matureHist);
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

    pNewDstrb->numPairMode = 0;
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

SR_Status SR_FragLenDstrbSetPairMode(SR_FragLenDstrb* pDstrb, const char* cmdArg)
{
    for (unsigned int i = 0; i != NUM_TOTAL_PAIR_MODE; ++i)
        pDstrb->validModeMap[i] = -1;

    if (cmdArg == NULL)
    {
        pDstrb->numPairMode = 2;
        pDstrb->validModeMap[1] = 0;
        pDstrb->validModeMap[5] = 0;

        return SR_OK;
    }

    unsigned int argLen = strlen(cmdArg);
    if (argLen != 7 || argLen != 3)
    {
        SR_ErrMsg("ERROR: Invalid pair mode argument.(should be either 2 or 4 unique numbers between 1-8, separated by comma)");
        return SR_ERR;
    }

    unsigned int pairModeSet = 0;
    for (unsigned int i = 0; i < argLen; i += 2)
    {
        int pairMode = cmdArg[i] - '0';
        ++(pDstrb->numPairMode);

        if (pairMode < 1 || pairMode > 8)
        {
            SR_ErrMsg("ERROR: Invalid pair mode: %c\n", cmdArg[i]);
            return SR_ERR;
        }
        else if (cmdArg[i + 1] != ',' && cmdArg[i + 1] != '\0')
        {
            SR_ErrMsg("ERROR: Invalid pair mode argument.(should be either 2 or 4 unique numbers between 1-8, separated by comma)");
            return SR_ERR;
        }
        else if (pDstrb->validModeMap[pairMode] != -1)
        {
            SR_ErrMsg("ERROR: Pair mode number should be unique. Found %d twice\n", pairMode);
            return SR_ERR;
        }

        if (pDstrb->numPairMode % 2 != 0)
        {
            pairModeSet = (pairMode << 3);
        }
        else
        {
            pairModeSet |= pairMode;
            if (SR_PairModeSetMap[pairModeSet] == 0)
            {
                SR_ErrMsg("ERROR: Pair modes %d and %d are not compatible.\n", pairModeSet >> 3, pairModeSet & 7);
                return SR_ERR;
            }
        }

        pDstrb->validModeMap[pairMode] = (pDstrb->numPairMode - 1) / 2;
    }

    return SR_OK;
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
        SR_ErrQuit("ERROR: Not enough memory for the storage of the histogram in the fragment length distribution object.\n");

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

    // if the pair mode is not valid
    // we only updated the count of the invalid pair and return
    int8_t histIndex = pDstrb->validModeMap[pPairStats->pairMode];
    if (histIndex < 0)
    {
        ++(pCurrHist->modeCount[INVALID_PAIR_MODE_SET_INDEX]);
        return SR_OK;
    }

    // because we can have up to 2 different pair mode sets (4 differen pair modes)
    // we should choose which histogram we should update
    // the first one (with index 0 or 1) or the second one(2, 3)
    khash_t(fragLen)* pCurrHash = pCurrHist->rawHist[histIndex];

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

    ++(pCurrHist->modeCount[histIndex]);
    pCurrHist->rawHist[histIndex] = pCurrHash;

    return SR_OK;
}

void SR_FragLenDstrbFinalize(SR_FragLenDstrb* pDstrb)
{
    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        SR_FragLenHistToMature(&(pDstrb->pHists[i]));
        SR_FragLenHistClearRaw(&(pDstrb->pHists[i]));
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
