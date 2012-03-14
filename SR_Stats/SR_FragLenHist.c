/*
 * =====================================================================================
 *
 *       Filename:  SR_FragLenHist.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/08/2012 02:24:40 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "khash.h"
#include "SR_Error.h"
#include "SR_Utilities.h"
#include "SR_FragLenHist.h"

// index of the mode count for invalid pair mode 
#define INVALID_PAIR_MODE_SET_INDEX 1

#define DEFAULT_NUM_HIST_ELMNT 200

// fragment length hash
KHASH_MAP_INIT_INT(fragLen, uint32_t);


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
    khash_t(fragLen)* pRawHist = pHist->rawHist;

    pHist->size = kh_size(pRawHist);

    if (pHist->size > pHist->capacity)
    {
        free(pHist->fragLen);
        free(pHist->freq);

        pHist->capacity = 2 * pHist->size;

        pHist->fragLen = (uint32_t*) malloc(pHist->capacity * sizeof(uint32_t));
        if(pHist->fragLen == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the fragment length array in the fragment length histogram object.\n");

        pHist->freq = (uint32_t*) malloc(pHist->capacity * sizeof(uint32_t));
        if(pHist->freq == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the frequency array in the fragment length histogram object.\n");
    }

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
    double cdf = 0;
    uint32_t fragLenQual = 0;

    SR_Bool foundMedian = FALSE;
    for (unsigned int j = 0; j != pHist->size; ++j)
    {
        khiter_t khIter = kh_get(fragLen, pRawHist, pHist->fragLen[j]);
        if (khIter == kh_end(pRawHist))
            SR_ErrQuit("ERROR: Cannot find the fragment length frequency from the hash table.\n");

        pHist->freq[j] = kh_value(pRawHist, khIter);

        totalFragLen += pHist->fragLen[j] * pHist->freq[j];
        cumFreq += pHist->freq[j];
        cdf = cumFreq / totalFreq;

        if (!foundMedian && cdf >= 0.5)
        {
            pHist->median = pHist->fragLen[j];
            foundMedian = TRUE;
        }

        cdf = cdf > 0.5 ? 1.0 - cdf : cdf;

        fragLenQual = DoubleRoundToInt(-10.0 * log10(cdf));
        kh_value(pRawHist, khIter) = fragLenQual;
    }

    pHist->mean = totalFragLen / totalFreq;

    pHist->stdev = 0.0;
    for (unsigned int j = 0; j != pHist->size; ++j)
    {
        pHist->stdev += (double) pHist->freq[j] * pow(pHist->mean - pHist->fragLen[j], 2);
    }

    if (totalFreq != 1)
        pHist->stdev = sqrt(pHist->stdev / (double) (totalFreq - 1));
}

SR_FragLenHistArray* SR_FragLenHistArrayAlloc(unsigned int capacity)
{
    SR_FragLenHistArray* pHistArray = NULL;
    SR_ARRAY_ALLOC(pHistArray, capacity, SR_FragLenHistArray, SR_FragLenHist);

    return pHistArray;
}

void SR_FragLenHistArrayFree(SR_FragLenHistArray* pHistArray)
{
    if (pHistArray != NULL)
    {
        for (unsigned int i = 0; i != pHistArray->capacity; ++i)
        {
            free(pHistArray->data[i].fragLen);
            free(pHistArray->data[i].freq);

            kh_destroy(fragLen, pHistArray->data[i].rawHist);
        }

        free(pHistArray->data);
        free(pHistArray);
    }
}

void SR_FragLenHistArrayClear(SR_FragLenHistArray* pHistArray)
{
    for (unsigned int i = 0; i != pHistArray->size; ++i)
    {
        pHistArray->data[i].size = 0;
        pHistArray->data[i].mean = 0.0;
        pHistArray->data[i].median = 0.0;
        pHistArray->data[i].stdev = 0.0;
        pHistArray->data[i].modeCount[0] = 0;
        pHistArray->data[i].modeCount[1] = 0;

        kh_clear(fragLen, pHistArray->data[i].rawHist);
    }
}

void SR_FragLenHistArrayInit(SR_FragLenHistArray* pHistArray, unsigned int newSize)
{
    SR_FragLenHistArrayClear(pHistArray);

    if (newSize > pHistArray->capacity)
    {
        SR_ARRAY_RESIZE(pHistArray, newSize * 2, SR_FragLenHist);
        memset(pHistArray->data + pHistArray->size, 0, (newSize * 2 - pHistArray->size) * sizeof(SR_FragLenHist));
    }

    for (unsigned int i = 0; i != newSize; ++i)
    {
        if (pHistArray->data[i].rawHist == NULL)
        {
            pHistArray->data[i].rawHist = kh_init(fragLen);
            kh_resize(fragLen, pHistArray->data[i].rawHist, DEFAULT_NUM_HIST_ELMNT);
        }
    }

    pHistArray->size = newSize;
}

SR_Status SR_FragLenHistArrayUpdate(SR_FragLenHistArray* pHistArray, unsigned int backHistIndex, uint32_t fragLen)
{
    if (backHistIndex > pHistArray->size)
        return SR_ERR;

    SR_FragLenHist* pCurrHist = pHistArray->data + (pHistArray->size - backHistIndex);

    // if the pair mode is not valid
    // we only updated the count of the invalid pair and return
    if (fragLen == 0)
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
    khiter_t khIter = kh_put(fragLen, pCurrHash, fragLen, &ret);

    if (ret == 0)
        kh_value(pCurrHash, khIter) += 1;
    else
        kh_value(pCurrHash, khIter) = 1;

    ++(pCurrHist->modeCount[0]);
    pCurrHist->rawHist = pCurrHash;

    return SR_OK;
}

int SR_FragLenHistArrayGetFragLenQual(const SR_FragLenHistArray* pHistArray, unsigned int backHistIndex, uint32_t fragLen)
{
    if (backHistIndex > pHistArray->size)
        return -1;

    SR_FragLenHist* pCurrHist = pHistArray->data + (pHistArray->size - backHistIndex);

    khash_t(fragLen)* pCurrHash = pCurrHist->rawHist;
    khiter_t khIter = kh_get(fragLen, pCurrHash, fragLen);

    int fragLenQual = -1;
    if (khIter != kh_end(pCurrHash))
        fragLenQual = kh_value(pCurrHash, khIter);

    return fragLenQual;
}

void SR_FragLenHistArrayFinalize(SR_FragLenHistArray* pHistArray)
{
    for (unsigned int i = 0; i != pHistArray->size; ++i)
        SR_FragLenHistToMature(&(pHistArray->data[i]));
}

void SR_FragLenHistArrayWrite(const SR_FragLenHistArray* pHistArray, FILE* output)
{
    for (unsigned int i = 0; i != pHistArray->size; ++i)
    {
        fwrite(&(pHistArray->data[i].size), sizeof(uint32_t), 1, output);
        fwrite(pHistArray->data[i].fragLen, sizeof(uint32_t), pHistArray->data[i].size, output);
        fwrite(pHistArray->data[i].freq, sizeof(uint32_t), pHistArray->data[i].size, output);
    }

    fflush(output);
}
