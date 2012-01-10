/*
 * =====================================================================================
 *
 *       Filename:  SR_ReadPairDetector.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/08/2011 11:45:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <math.h>

#include "SR_Utilities.h"
#include "SR_Error.h"
#include "SR_ReadPairDetector.h"

#define SR_DEFAULT_RP_CAPACITY 50

static const SV_EventType SV_EventTypeMap[8][8] = 
{
    {SV_UNKNOWN, SV_INVERSION3, SV_INVERSION5, SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN},
    {SV_INVERSION3, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION5, SV_UNKNOWN, SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN},
    {SV_INVERSION5, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION3, SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN, SV_UNKNOWN},
    {SV_UNKNOWN, SV_INVERSION5, SV_INVERSION3, SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN, SV_TANDEM_DUP},
    {SV_TANDEM_DUP, SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION3, SV_INVERSION5, SV_UNKNOWN},
    {SV_UNKNOWN, SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN, SV_INVERSION3, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION5},
    {SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION5, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION3},
    {SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN, SV_INVERSION5, SV_INVERSION3, SV_UNKNOWN}
};

static void SR_ReadPairLoadProb(SR_ReadPairInfo* pInfo, const SR_FragLenDstrb* pDstrb, unsigned int begin, unsigned int end)
{
    const uint32_t* fragLenArray = pDstrb->pHists[pInfo->readGrpID].fragLen;
    const double* probArray = pDstrb->pHists[pInfo->readGrpID].cdf;

    unsigned int min = begin;
    unsigned int max = end;


    while (min <= max)
    {
        unsigned int mid = (min + max) / 2;

        if (pInfo->fragLen == fragLenArray[mid])
        {
            double prob = begin == 0 ? probArray[mid] : (1 - probArray[mid]);
            if (prob == 0.0)
            {
                pInfo->probPower = 255;
                pInfo->probNum = 1;
            }
            else
            {
                pInfo->probPower = -ceil(log10(prob));
                pInfo->probNum = (int32_t) prob * pow(10, pInfo->probPower);
            }

            return;
        }
        else if (pInfo->fragLen > fragLenArray[mid])
        {
            min = mid + 1;
        }
        else
        {
            if (mid == 0)
                break;

            max = mid - 1;
        }
    }

    pInfo->probPower = 255;
    pInfo->probNum = 1;
}

static int CompareAttrbt(const void* a, const void* b)
{
    const SR_ReadPairAttrbt* first = a;
    const SR_ReadPairAttrbt* second = b;

    if (first->firstAttribute < second->firstAttribute)
    {
        return -1;
    }
    else if (first->firstAttribute > second->firstAttribute)
    {
        return 1;
    }
    else
    {
        if (first->secondAttribute < second->secondAttribute)
        {
            return -1;
        }
        else if (first->secondAttribute > second->secondAttribute)
        {
            return 1;
        }
        else
        {
            if (first->readGrpID <= second->readGrpID)
            {
                return -1;
            }
            else
            {
                return 1;
            }
        }
    }
}

SR_ReadPairInfoTable* SR_ReadPairInfoTableAlloc(uint32_t numChr, uint32_t numRG, double cutoff)
{
    SR_ReadPairInfoTable* pNewInfoTable = (SR_ReadPairInfoTable*) malloc(sizeof(SR_ReadPairInfoTable));
    if (pNewInfoTable == NULL)
        SR_ErrQuit("ERROR: Not enough memory for read pair information table.\n");

    for (unsigned int i = 0; i != SR_NUM_SV_TYPES; ++i)
    {
        SR_ARRAY_INIT(&(pNewInfoTable->arrays[i]), SR_DEFAULT_RP_CAPACITY, SR_ReadPairInfo);
    }

    pNewInfoTable->numChr = numChr;
    pNewInfoTable->chrIndex = (uint32_t*) calloc(SR_NUM_SV_TYPES * numChr, sizeof(uint32_t));
    if (pNewInfoTable->chrIndex == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of chromosome index in the read pair information table.\n");

    pNewInfoTable->numRG = numRG;

    return pNewInfoTable;
}

void SR_ReadPairInfoTableFree(SR_ReadPairInfoTable* pInfoTable)
{
    if (pInfoTable != NULL)
    {
        free(pInfoTable->chrIndex);
        free(pInfoTable);
    }
}

SR_ReadPairAttrbtArray* SR_ReadPairAttrbtArrayAlloc(uint32_t dataCap, uint32_t boundCap)
{
    SR_ReadPairAttrbtArray* pAttrbtArray = (SR_ReadPairAttrbtArray*) malloc(sizeof(SR_ReadPairAttrbtArray));
    if (pAttrbtArray == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a read pair attribute array object.\n");

    pAttrbtArray->data = (SR_ReadPairAttrbt*) malloc(dataCap * sizeof(SR_ReadPairAttrbt));
    if (pAttrbtArray->data == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the read pair attributes in the read pair attribute array object.\n");

    pAttrbtArray->dataSize = 0;
    pAttrbtArray->dataCap = dataCap;

    pAttrbtArray->bound = (double*) malloc(boundCap * sizeof(double));
    if (pAttrbtArray->bound == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the cutoff in the read pair attribute array object.\n");

    pAttrbtArray->boundSize = 0;
    pAttrbtArray->boundCap = boundCap;

    return pAttrbtArray;
}

void SR_ReadPairAttrbtArrayFree(SR_ReadPairAttrbtArray* pAttrbtArray)
{
    if (pAttrbtArray != NULL)
    {
        free(pAttrbtArray->data);
        free(pAttrbtArray->bound);
        
        free(pAttrbtArray);
    }
}

SR_Status SR_ReadPairInfoTableUpdate(SR_ReadPairInfoTable* pInfoTable, const SR_BamNode* pUpAlgn, const SR_BamNode* pDownAlgn, const SR_FragLenDstrb* pDstrb)
{
    SR_ReadPairInfo info;

    info.upPos = pUpAlgn->alignment.core.pos;
    info.downPos = pDownAlgn->alignment.core.pos;

    // get the statistics of the read pair
    SR_BamPairStats pairStats;
    SR_LoadPairStats(&pairStats, pUpAlgn);

    info.pairMode = pairStats.pairMode;
    SR_FragLenDstrbGetRGIndex(&(info.readGrpID), pDstrb, pairStats.RG);

    info.eventType = SV_UNKNOWN;
    info.probNum = 0;
    info.probPower = 0;

    info.downRefID = pDownAlgn->alignment.core.tid;
    info.upRefID = pUpAlgn->alignment.core.tid;

    if (pUpAlgn->alignment.core.tid == pDownAlgn->alignment.core.tid)
        info.fragLen = abs(pUpAlgn->alignment.core.isize);
    else
    {
        info.fragLen = -1;
        info.eventType = SV_INTER_CHR_TRNSLCTN;
        return SR_OK;
    }

    if (SR_IsValidPairMode(pDstrb, pairStats.pairMode))
    {
        // any reads with fragment length at the edge of the fragment length distribution will be kept as SV candidates
        uint32_t lowerCutoffIndex = SR_GetHistCutoffIndex(pDstrb, info.readGrpID, DSTRB_LOWER_CUTOFF);
        uint32_t upperCutoffIndex = SR_GetHistCutoffIndex(pDstrb, info.readGrpID, DSTRB_UPPER_CUTOFF);

        uint32_t lowerCutoff = SR_GetHistCutoffValue(pDstrb, info.readGrpID, lowerCutoffIndex);
        uint32_t upperCutoff = SR_GetHistCutoffValue(pDstrb, info.readGrpID, upperCutoffIndex);

        if (pairStats.fragLen < lowerCutoff)
        {
            SR_ReadPairLoadProb(&info, pDstrb, 0, lowerCutoffIndex);
            info.eventType = SV_INSERTION;
        }
        else if (pairStats.fragLen > upperCutoff)
        {
            SR_ReadPairLoadProb(&info, pDstrb, upperCutoffIndex, pDstrb->pHists[info.readGrpID].size - 1);
            info.eventType = SV_DELETION;
        }
    }
    else
    {
        SV_EventType eventFirst = SV_EventTypeMap[pairStats.pairMode][pDstrb->validMode[0]];
        SV_EventType eventSecond = SV_EventTypeMap[pairStats.pairMode][pDstrb->validMode[1]];

        if ((eventFirst == SV_UNKNOWN && eventSecond == SV_UNKNOWN) || (eventFirst != SV_UNKNOWN && eventSecond != SV_UNKNOWN))
        {
            info.eventType = SV_UNKNOWN;
        }
        else if (eventFirst != SV_UNKNOWN)
        {
            info.eventType = eventFirst;
        }
        else
        {
            info.eventType = eventSecond;
        }
    }

    if (info.eventType != SV_UNKNOWN)
    {
        SR_ReadPairInfoArray* pInfoArray = &(pInfoTable->arrays[info.eventType - 1]);
        SR_ARRAY_PUSH(pInfoArray, &info, SR_ReadPairInfo);

        unsigned int index = (info.eventType - 1) * pInfoTable->numChr + info.downRefID;
        ++(pInfoTable->chrIndex[index]);
    }
    else
        return SR_NOT_FOUND;

    return SR_OK;
}

void SR_ReadPairAttrbtArrayResize(SR_ReadPairAttrbtArray* pAttrbtArray, uint32_t newDataCap, uint32_t newBoundCap)
{
    pAttrbtArray->dataSize = 0;
    pAttrbtArray->boundSize = 0;

    if (newDataCap > pAttrbtArray->dataCap)
    {
        free(pAttrbtArray->data);
        pAttrbtArray->data = (SR_ReadPairAttrbt*) malloc(newDataCap * sizeof(SR_ReadPairAttrbt));
        if (pAttrbtArray->data == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the read pair attributes in the read pair attribute array object.\n");

        pAttrbtArray->dataCap = newDataCap;
    }

    if (newBoundCap > pAttrbtArray->boundCap)
    {
        free(pAttrbtArray->bound);
        pAttrbtArray->bound = (double*) malloc(newDataCap * sizeof(double));
        if (pAttrbtArray->bound == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the boundaries in the read pair attribute array object.\n");

        pAttrbtArray->boundCap = newDataCap;
    }
}

void SR_ReadPairMake(SR_ReadPairAttrbtArray* pAttrbtArray, const SR_FragLenDstrb* pDstrb, const SR_ReadPairInfoArray* pInfoArray, SV_EventType eventType)
{
    if (pInfoArray->size == 0)
        return;

    SR_ReadPairAttrbtArrayResize(pAttrbtArray, pInfoArray->size, pDstrb->size * 2);

    pAttrbtArray->eventType = eventType;
    pAttrbtArray->dataSize = pInfoArray->size;
    pAttrbtArray->boundSize = pDstrb->size * 2;

    if (eventType == SV_INTER_CHR_TRNSLCTN)
    {
        for (unsigned int i = 0; i != pInfoArray->size; ++i)
        {
            pAttrbtArray->data[i].origIdex = i;
            pAttrbtArray->data[i].readGrpID = pInfoArray->data[i].readGrpID;

            pAttrbtArray->data[i].firstAttribute = pInfoArray->data[i].upPos;
            pAttrbtArray->data[i].secondAttribute = pInfoArray->data[i].downRefID * 1e10 + pInfoArray->data[i].downPos;
        }

        for (unsigned int i = 0; i != pAttrbtArray->boundSize; i += 2)
        {
            uint32_t lowerCutoffIndex = SR_GetHistCutoffIndex(pDstrb, i / 2, DSTRB_LOWER_CUTOFF);
            uint32_t upperCutoffIndex = SR_GetHistCutoffIndex(pDstrb, i / 2, DSTRB_UPPER_CUTOFF);

            uint32_t lowerCutoff = SR_GetHistCutoffValue(pDstrb, i / 2, lowerCutoffIndex);
            uint32_t upperCutoff = SR_GetHistCutoffValue(pDstrb, i / 2, upperCutoffIndex);

            pAttrbtArray->bound[i] = upperCutoff - lowerCutoff;
            pAttrbtArray->bound[i + 1] = upperCutoff - lowerCutoff;
        }
    }
    else if (eventType <= SV_TANDEM_DUP)
    {

        for (unsigned int i = 0; i != pInfoArray->size; ++i)
        {
            pAttrbtArray->data[i].origIdex = i;
            pAttrbtArray->data[i].readGrpID = pInfoArray->data[i].readGrpID;

            double median = SR_GetHistMedian(pDstrb, pInfoArray->data[i].readGrpID);

            pAttrbtArray->data[i].firstAttribute = pInfoArray->data[i].upPos + (double) pInfoArray->data[i].fragLen / 2;
            pAttrbtArray->data[i].secondAttribute = pInfoArray->data[i].fragLen - median;
        }

        for (unsigned int i = 0; i != pAttrbtArray->boundSize; i += 2)
        {
            uint32_t lowerCutoffIndex = SR_GetHistCutoffIndex(pDstrb, i / 2, DSTRB_LOWER_CUTOFF);
            uint32_t upperCutoffIndex = SR_GetHistCutoffIndex(pDstrb, i / 2, DSTRB_UPPER_CUTOFF);

            uint32_t lowerCutoff = SR_GetHistCutoffValue(pDstrb, i / 2, lowerCutoffIndex);
            uint32_t upperCutoff = SR_GetHistCutoffValue(pDstrb, i / 2, upperCutoffIndex);

            double median = SR_GetHistMedian(pDstrb, i / 2);

            pAttrbtArray->bound[i] = median;
            pAttrbtArray->bound[i + 1] = upperCutoff - lowerCutoff;

            if (eventType == SV_INVERSION3 || eventType == SV_INVERSION5)
                pAttrbtArray->bound[i + 1] = 2 * median;
        }
    }
    else
    {
        //TODO: MEI
    }

    qsort(pAttrbtArray->data, pAttrbtArray->dataSize, sizeof(pAttrbtArray->data[0]), CompareAttrbt);
}
