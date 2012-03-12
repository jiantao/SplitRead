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
#include "SR_Utilities.h"
#include "SR_FragLenDstrb.h"


//===============================
// Type and constant definition
//===============================

// default capacity of fragment length distribution
#define DEFAULT_RG_CAP 10

#define DEFAULT_SM_CAP 20

// index of the mode count for invalid pair mode 
#define INVALID_PAIR_MODE_SET_INDEX 1

// read group name hash
KHASH_MAP_INIT_STR(readGrp, uint32_t);

// sample name hash
KHASH_MAP_INIT_STR(sample, uint32_t);

// fragment length hash
KHASH_MAP_INIT_INT(fragLen, uint32_t);

//===================
// Static functions
//===================

static SR_Status SR_FragLenDstrbAddSample(int* pSampleID, SR_FragLenDstrb* pDstrb, const char* sampleNamePos)
{
    const char* sampleNameEnd = strpbrk(sampleNamePos, " \t\n\0");
    int sampleNameLen = sampleNameEnd - sampleNamePos;
    if (sampleNameLen > 0)
    {
        char* buff = (char*) malloc((sampleNameLen + 1) * sizeof(char));
        if (buff == NULL)
            SR_ErrQuit("ERROR: Not enought memory for the read group ID in the fragment length distribution object.\n");

        buff[sampleNameLen] = '\0';
        memcpy(buff, sampleNamePos, sampleNameLen);

        int ret = 0;

        khash_t(sample)* pSampleHash = pDstrb->pSampleHash;
        khiter_t khIter = kh_put(sample, pSampleHash, buff, &ret);

        if (ret == 0)
        {
            free(buff);
            *pSampleID = kh_value(pSampleHash, khIter);
        }
        else
        {
            if (pDstrb->sizeSM == pDstrb->capacitySM)
            {
                pDstrb->capacitySM *= 2;
                pDstrb->pSamples = (char**) realloc(pDstrb->pSamples, pDstrb->capacitySM * sizeof(char*));
                if (pDstrb->pSamples == NULL)
                    SR_ErrQuit("ERROR: Not enought memory for the sample names in the fragment length distribution object.\n");
            }

            *pSampleID = pDstrb->sizeSM;
            kh_value(pSampleHash, khIter) = pDstrb->sizeSM;
            pDstrb->pSamples[pDstrb->sizeSM] = buff;
            ++(pDstrb->sizeSM);
        }

        return SR_OK;
    }

    return SR_ERR;
}

static SR_Status SR_FragLenDstrbAddReadGrp(SR_FragLenDstrb* pDstrb, const char* readGrpNamePos, int sampleID)
{
    // expand the array if necessary
    if (pDstrb->sizeRG == pDstrb->capacityRG)
    {
        pDstrb->capacityRG *= 2;
        pDstrb->pReadGrps = (char**) realloc(pDstrb->pReadGrps, pDstrb->capacityRG * sizeof(char*));
        if (pDstrb->pReadGrps == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of read group names in the fragment length distribution object.\n");

        pDstrb->pSampleMap = (int32_t*) realloc(pDstrb->pSampleMap, pDstrb->capacityRG * sizeof(int32_t));
        if (pDstrb->pSampleMap == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the sample ID map in the fragment length distribution object.\n");
    }

    const char* readGrpNameEnd = strpbrk(readGrpNamePos, " \t\n\0");
    size_t readGrpNameLen = readGrpNameEnd - readGrpNamePos;
    if (readGrpNameLen > 0)
    {
        pDstrb->pReadGrps[pDstrb->sizeRG] = (char*) calloc(readGrpNameLen + 1, sizeof(char));
        if (pDstrb->pReadGrps[pDstrb->sizeRG] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");

        memcpy(pDstrb->pReadGrps[pDstrb->sizeRG], readGrpNamePos, readGrpNameLen);

        int ret = 0;

        khash_t(readGrp)* pReadGrpHash = pDstrb->pReadGrpHash;
        khiter_t khIter = kh_put(readGrp, pReadGrpHash, pDstrb->pReadGrps[pDstrb->sizeRG], &ret);

        if (ret != 0)
        {
            pDstrb->pSampleMap[pDstrb->sizeRG] = sampleID;
            kh_value(pReadGrpHash, khIter) = pDstrb->sizeRG;
            ++(pDstrb->sizeRG);

            return SR_OK;
        }
        else
        {
            free(pDstrb->pReadGrps[pDstrb->sizeRG]);
            pDstrb->pReadGrps[pDstrb->sizeRG] = NULL;
            SR_ErrMsg("ERROR: Found a duplicated read group ID.\n");
        }
    }

    return SR_ERR;
}

static inline void SR_FragLenHistClearRaw(SR_FragLenHist* pHist)
{
    kh_destroy(fragLen, pHist->rawHist);
    pHist->rawHist = NULL;
}

static inline void SR_FragLenHistClear(SR_FragLenHist* pHist)
{
    SR_FragLenHistClearRaw(pHist);

    free(pHist->fragLen);
    free(pHist->freq);

    pHist->fragLen = NULL;
    pHist->freq = NULL;
}

static inline void SR_FragLenHistClearHist(SR_FragLenDstrb* pDstrb)
{
    for (unsigned int i = 0; i != pDstrb->sizeHist; ++i)
        SR_FragLenHistClear(&(pDstrb->pHists[i]));

    pDstrb->sizeHist = 0;
}

static inline int CompareFragLenBin(const void* a, const void* b)
{
    const uint32_t* first = a;
    const uint32_t* second = b;

    if (*first < *second)
        return -1;
    else if (*first > *second)
        return 1;
    else
        return 0;
}

static void SR_FragLenHistToMature(SR_FragLenHist* pHist)
{
    const khash_t(fragLen)* pRawHist = pHist->rawHist;

    pHist->size = kh_size(pRawHist);
    pHist->fragLen = (uint32_t*) malloc(pHist->size * sizeof(uint32_t));
    if(pHist->fragLen == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the fragment length array in the fragment length histogram object.\n");

    pHist->freq = (uint32_t*) malloc(pHist->size * sizeof(uint32_t));
    if(pHist->freq == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the frequency array in the fragment length histogram object.\n");

    unsigned int i = 0;
    for (khiter_t khIter = kh_begin(pRawHist); khIter != kh_end(pRawHist); ++khIter)
    {
        if (kh_exist(pRawHist, khIter))
        {
            pHist->fragLen[i] = kh_key(pRawHist, khIter);
            ++i;
        }
    }

    qsort(pHist->fragLen, pHist->size, sizeof(uint32_t), CompareFragLenBin);

    double cumFreq = 0.0;
    double totalFragLen = 0.0;
    uint64_t totalFreq = pHist->modeCount[0];

    for (unsigned int j = 0; j != pHist->size; ++j)
    {
        khiter_t khIter = kh_get(fragLen, pHist->rawHist, pHist->fragLen[j]);
        if (khIter == kh_end((khash_t(fragLen)*) pHist->rawHist))
            SR_ErrQuit("ERROR: Cannot find the fragment length frequency from the hash table.\n");

        pHist->freq[j] = kh_value((khash_t(fragLen)*) pHist->rawHist, khIter);

        totalFragLen += pHist->fragLen[j] * pHist->freq[j];
        cumFreq += pHist->freq[j];
    }

    pHist->mean = totalFragLen / totalFreq;
    pHist->median = FindMedianUint(pHist->fragLen, pHist->size);

    pHist->stdev = 0.0;
    for (unsigned int j = 0; j != pHist->size; ++j)
    {
        pHist->stdev += (double) pHist->freq[j] * pow(pHist->mean - pHist->fragLen[j], 2);
    }

    if (totalFreq != 1)
        pHist->stdev = sqrt(pHist->stdev / (double) (totalFreq - 1));
}


//===============================
// Constructors and Destructors
//===============================

SR_FragLenDstrb* SR_FragLenDstrbAlloc(void)

{
    SR_FragLenDstrb* pNewDstrb = (SR_FragLenDstrb*) calloc(1, sizeof(*pNewDstrb));
    if (pNewDstrb == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the fragment length distribution object.\n");

    pNewDstrb->pSamples = (char**) malloc(sizeof(char*) * DEFAULT_SM_CAP);
    if (pNewDstrb->pSamples == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the sample names in the fragment length distribution object.\n");

    pNewDstrb->pReadGrps = (char**) malloc(sizeof(char*) * DEFAULT_RG_CAP);
    if (pNewDstrb->pReadGrps == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the read group names in the fragment length distribution object.\n");

    pNewDstrb->pSeqTech = (int8_t*) malloc(sizeof(int8_t) * DEFAULT_RG_CAP);
    if (pNewDstrb->pReadGrps == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the sequencing technology in the fragment length distribution object.\n");

    pNewDstrb->pSampleMap = (int32_t*) malloc(sizeof(int32_t) * DEFAULT_RG_CAP);
    if (pNewDstrb->pSampleMap == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the sample-read-group-ID map in the fragment length distribution object.\n");

    pNewDstrb->pReadGrpHash = kh_init(readGrp);
    kh_resize(readGrp, pNewDstrb->pReadGrpHash, DEFAULT_RG_CAP);

    pNewDstrb->pSampleHash = kh_init(sample);
    kh_resize(sample, pNewDstrb->pSampleHash, DEFAULT_SM_CAP);

    pNewDstrb->pHists = NULL;

    pNewDstrb->sizeRG = 0;
    pNewDstrb->capacityRG = DEFAULT_RG_CAP;
    pNewDstrb->sizeSM = 0;
    pNewDstrb->capacitySM = DEFAULT_SM_CAP;
    pNewDstrb->sizeHist = 0;
    pNewDstrb->capacityHist = 0;

    return pNewDstrb;
}

void SR_FragLenDstrbFree(SR_FragLenDstrb* pDstrb)
{
    if (pDstrb != NULL)
    {
        if (pDstrb->pHists != NULL)
        {
            for (unsigned int i = 0; i != pDstrb->sizeHist; ++i) 
                SR_FragLenHistClear(&(pDstrb->pHists[i]));

            free(pDstrb->pHists);
        }

        if (pDstrb->pReadGrps != NULL)
        {
            for (unsigned int i = 0; i != pDstrb->sizeRG; ++i)
                free(pDstrb->pReadGrps[i]);

            free(pDstrb->pReadGrps);
        }

        if (pDstrb->pSamples != NULL)
        {
            for (unsigned int i = 0; i != pDstrb->sizeSM; ++i)
                free(pDstrb->pSamples[i]);

            free(pDstrb->pSamples);
        }

        free(pDstrb->pSampleMap);

        kh_destroy(readGrp, pDstrb->pReadGrpHash);
        kh_destroy(sample, pDstrb->pSampleHash);

        free(pDstrb);
    }
}


//======================
// Interface functions
//======================

SR_PairMode SR_GetPairMode(const bam1_t* pAlignment)
{
    unsigned int upMode = 0;
    unsigned int  downMode = 0;

    if ((pAlignment->core.flag & BAM_FREVERSE) != 0)
    {
        upMode |= 1;
    }

    if ((pAlignment->core.flag & BAM_FMREVERSE) != 0)
    {
        downMode |= 1;
    }

    if ((pAlignment->core.flag & BAM_FREAD1) != 0 && (pAlignment->core.flag & BAM_FREAD2) != 0)
        return SR_BAD_PAIR_MODE;
    if ((pAlignment->core.flag & BAM_FREAD1) != 0)
        downMode |= (1 << 1);
    else if ((pAlignment->core.flag & BAM_FREAD2))
        upMode |= (1 << 1);
    else
        return SR_BAD_PAIR_MODE;

    if (pAlignment->core.isize < 0)
        SR_SWAP(upMode, downMode, unsigned int);

    return ((SR_PairMode) SR_PairModeMap[(upMode << 2) | downMode]);
}

SR_Status SR_LoadPairStats(SR_BamPairStats* pPairStats, const bam1_t* pAlignment)
{
    pPairStats->pairMode = SR_GetPairMode(pAlignment);
    if (pPairStats->pairMode == SR_BAD_PAIR_MODE)
        return SR_ERR;

    pPairStats->fragLen = abs(pAlignment->core.isize);

    static const char tagRG[2] = {'R', 'G'};
    uint8_t* rgPos = bam_aux_get(pAlignment, tagRG);
    if (rgPos != NULL)
        pPairStats->RG = bam_aux2Z(rgPos);
    else
        pPairStats->RG = NULL;

    return SR_OK;
}

void SR_FragLenDstrbSetPairMode(SR_FragLenDstrb* pDstrb, const int8_t* pValidPairMode)
{
    pDstrb->validMode[0] = pValidPairMode[0];
    pDstrb->validMode[1] = pValidPairMode[1];

    pDstrb->validModeMap &= 0;
    pDstrb->validModeMap |= (1 << (pValidPairMode[0]));
    pDstrb->validModeMap |= (1 << (pValidPairMode[1]));
}

SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader)
{
    SR_Status status = SR_OK;

    const char* tagPos = pBamHeader->pOrigHeader->text;
    while ((tagPos = strstr(tagPos, "@RG")) != NULL)
    {
        const char* lineEnd = strpbrk(tagPos, "\n");
        const char* sampleNamePos = strstr(tagPos, "SM:");
        int32_t sampleID = 0;

        if (sampleNamePos != NULL && sampleNamePos < lineEnd)
        {
            sampleNamePos += 3;
            status = SR_FragLenDstrbAddSample(&sampleID, pDstrb, sampleNamePos);
        }

        if (sampleNamePos == NULL || sampleNamePos >= lineEnd || status != SR_OK)
        {
            SR_ErrMsg("ERROR: the \"SM\" field is not found under the read group tag in the bam header.\n");

            ++tagPos;
            continue;
        }

        // get the name of the current read group
        const char* readGrpNamePos = strstr(tagPos, "ID:");
        if (readGrpNamePos != NULL && readGrpNamePos < lineEnd)
        {
            readGrpNamePos += 3;
            status = SR_FragLenDstrbAddReadGrp(pDstrb, readGrpNamePos, sampleID);
        }

        if (readGrpNamePos == NULL || readGrpNamePos >= lineEnd || status != SR_OK)
        {
            SR_ErrMsg("ERROR: the \"ID\" field is required under the read group tag.\n");

            ++tagPos;
            continue;
        }

        ++(pDstrb->sizeHist);
        ++tagPos;
    }

    if (pDstrb->sizeHist == 0)
    {
        SR_ErrMsg("ERROR: Read group tag is not found.\n");
        return SR_ERR;
    }

    // allocate the memory for the histograms
    if (pDstrb->sizeHist > pDstrb->capacityHist)
    {
        free(pDstrb->pHists);

        if (pDstrb->capacityHist != 0)
            pDstrb->capacityHist = 2 * pDstrb->sizeHist;
        else
            pDstrb->capacityHist = pDstrb->sizeHist;

        pDstrb->pHists = (SR_FragLenHist*) calloc(pDstrb->capacityHist, sizeof(SR_FragLenHist));
        if (pDstrb->pHists == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the histogram in the fragment length distribution object.\n");
    }

    return SR_OK;
}

SR_Status SR_FragLenDstrbGetRGIndex(int32_t* pReadGrpIndex, const SR_FragLenDstrb* pDstrb, const char* pReadGrpName)
{
    *pReadGrpIndex = 0;

    khash_t(readGrp)* pRgHash = pDstrb->pReadGrpHash;
    khiter_t khIter = kh_get(readGrp, pRgHash, pReadGrpName);
    if (khIter != kh_end(pRgHash))
    {
        *pReadGrpIndex = kh_value(pRgHash, khIter);
    }
    else
    {
        SR_ErrMsg("ERROR: Found a read group name that is not record in the header.\n");
        return SR_ERR;
    }

    return SR_OK;
}

SR_Status SR_FragLenDstrbUpdate(SR_FragLenDstrb* pDstrb, const SR_BamPairStats* pPairStats)
{
    int32_t dstrbIndex = 0;
    SR_Status RGStatus = SR_FragLenDstrbGetRGIndex(&dstrbIndex, pDstrb, pPairStats->RG);
    if (RGStatus == SR_ERR)
        return SR_ERR;

    dstrbIndex = pDstrb->sizeRG - dstrbIndex;
    // duplicated read group, probably because of redundant files
    if (dstrbIndex > pDstrb->sizeHist)
        return SR_ERR;

    dstrbIndex = pDstrb->sizeHist - dstrbIndex;
    SR_FragLenHist* pCurrHist = pDstrb->pHists + dstrbIndex;

    // if the pair mode is not valid
    // we only updated the count of the invalid pair and return
    if (SR_IsValidPairMode(pDstrb, pPairStats->pairMode))
    {
        ++(pCurrHist->modeCount[INVALID_PAIR_MODE_SET_INDEX]);
        return SR_OK;
    }

    // because we can have up to 2 different pair mode sets (4 differen pair modes)
    // we should choose which histogram we should update
    // the first one (with index 0 or 1) or the second one(2, 3)
    khash_t(fragLen)* pCurrHash = pCurrHist->rawHist;

    if (pCurrHash == NULL)
    {
        pCurrHash = kh_init(fragLen);
        kh_resize(fragLen, pCurrHash, 20);
    }

    int ret = 0;
    khiter_t khIter = kh_put(fragLen, pCurrHash, pPairStats->fragLen, &ret);

    if (ret == 0)
        kh_value(pCurrHash, khIter) += 1;
    else
        kh_value(pCurrHash, khIter) = 1;

    ++(pCurrHist->modeCount[0]);
    pCurrHist->rawHist = pCurrHash;

    return SR_OK;
}

void SR_FragLenDstrbFinalize(SR_FragLenDstrb* pDstrb)
{
    for (unsigned int i = 0; i != pDstrb->sizeHist; ++i)
    {
        SR_FragLenHistToMature(&(pDstrb->pHists[i]));
        SR_FragLenHistClearRaw(&(pDstrb->pHists[i]));
    }
}

void SR_FragLenDstrbInitHistHeader(const SR_FragLenDstrb* pDstrb, FILE* histOutput)
{
    uint32_t numHists = 0;
    unsigned int writeSize = fwrite(&numHists, sizeof(uint32_t), 1, histOutput);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the total number of histograms into the histogram file.\n");

    fflush(histOutput);
}

void SR_FragLenDstrbWriteHistHeader(uint32_t numHists, FILE* histOutput)
{
    fseeko(histOutput, 0, SEEK_SET);

    unsigned int writeSize = fwrite(&numHists, sizeof(uint32_t), 1, histOutput);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the total number of histograms into the histogram file.\n");

    fflush(histOutput);
}

void SR_FragLenDstrbWriteHist(const SR_FragLenDstrb* pDstrb, FILE* histOutput)
{
    unsigned int writeSize = 0;
    for (unsigned int i = 0; i != pDstrb->sizeHist; ++i)
    {
        writeSize = fwrite(&(pDstrb->pHists[i].size), sizeof(uint32_t), 1, histOutput);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Cannot write the size of the histogram into the histogram file.\n");

        writeSize = fwrite(pDstrb->pHists[i].fragLen, sizeof(uint32_t), pDstrb->pHists[i].size, histOutput);
        if (writeSize != pDstrb->pHists[i].size)
            SR_ErrQuit("ERROR: Cannot write the fragment length into the histogram file.\n");

        writeSize = fwrite(pDstrb->pHists[i].freq, sizeof(uint32_t), pDstrb->pHists[i].size, histOutput);
        if (writeSize != pDstrb->pHists[i].size)
            SR_ErrQuit("ERROR: Cannot write the fragment length frequency into the histogram file.\n");
    }

    fflush(histOutput);
}

void SR_FragLenDstrbWriteInfo(const SR_FragLenDstrb* pDstrb, FILE* libFile)
{
    unsigned int writeSize = 0;

    writeSize = fwrite(&(pDstrb->sizeSM), sizeof(uint32_t), 1, libFile);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the number of samples into the information file.\n");

    writeSize = fwrite(&(pDstrb->sizeRG), sizeof(uint32_t), 1, libFile);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the number of read groups into the information file.\n");

    for (unsigned int i = 0; i != pDstrb->sizeSM; ++i)
    {
        uint32_t sampleNameLen = strlen(pDstrb->pSamples[i]);
        writeSize = fwrite(&(sampleNameLen), sizeof(uint32_t), 1, libFile);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Cannot write the sample name length into the information file.\n");

        writeSize = fwrite(pDstrb->pSamples[i], sizeof(char), sampleNameLen, libFile);
        if (writeSize != sampleNameLen)
            SR_ErrQuit("ERROR: Cannot write the sample name into the information file.\n");
    }

    writeSize = fwrite(pDstrb->pSampleMap, sizeof(int32_t), pDstrb->sizeRG, libFile);
    if (writeSize != pDstrb->sizeRG)
        SR_ErrQuit("ERROR: Cannot write the read-group-to-sample map into the information file.\n");

    for (unsigned int i = 0; i != pDstrb->sizeSM; ++i)
    {
        uint32_t readGrpNameLen = strlen(pDstrb->pReadGrps[i]);
        writeSize = fwrite(&(readGrpNameLen), sizeof(uint32_t), 1, libFile);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Cannot write the readGrp name length into the information file.\n");

        writeSize = fwrite(pDstrb->pReadGrps[i], sizeof(char), readGrpNameLen, libFile);
        if (writeSize != readGrpNameLen)
            SR_ErrQuit("ERROR: Cannot write the readGrp name into the information file.\n");
    }

    fflush(libFile);
}
