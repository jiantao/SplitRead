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
    {SV_UNKNOWN, SV_INVERSION, SV_INVERSION, SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN},
    {SV_INVERSION, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION, SV_UNKNOWN, SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN},
    {SV_INVERSION, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION, SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN, SV_UNKNOWN},
    {SV_UNKNOWN, SV_INVERSION, SV_INVERSION, SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN, SV_TANDEM_DUP},
    {SV_TANDEM_DUP, SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION, SV_INVERSION, SV_UNKNOWN},
    {SV_UNKNOWN, SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN, SV_INVERSION, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION},
    {SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION, SV_UNKNOWN, SV_UNKNOWN, SV_INVERSION},
    {SV_UNKNOWN, SV_UNKNOWN, SV_UNKNOWN, SV_TANDEM_DUP, SV_UNKNOWN, SV_INVERSION, SV_INVERSION, SV_UNKNOWN}
};

static void SR_ReadPairLoadProb(SR_ReadPairInfo* pInfo, const SR_FragLenDstrb* pDstrb, unsigned int begin, unsigned int end)
{
    unsigned int histIndex = pDstrb->validModeMap[pInfo->pairMode];
    const uint32_t* fragLenArray = pDstrb->pHists[pInfo->readGrpID].fragLen[histIndex];
    const double* probArray = pDstrb->pHists[pInfo->readGrpID].cdf[histIndex];

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

SR_ReadPairInfoTable* SR_ReadPairInfoTableAlloc(uint32_t numChr, uint32_t numRG, double detectCutoff, double clusterCutoff)
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
    pNewInfoTable->detectBound = (uint32_t*) calloc(numRG * NUM_ALLOWED_HIST * 2, sizeof(uint32_t));
    if (pNewInfoTable->detectBound == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of detection boundaries in the read pair information table.\n");

    pNewInfoTable->clusterBound = (uint32_t*) calloc(numRG * NUM_ALLOWED_HIST * 2, sizeof(uint32_t));
    if (pNewInfoTable->clusterBound == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of cluster range in the read pair information table.\n");

    if (detectCutoff >= clusterCutoff)
        SR_ErrMsg("WARNING: The probablity of detection cutoff is too large(%f).\n", detectCutoff);

    pNewInfoTable->detectCutoff = detectCutoff;
    pNewInfoTable->clusterCutoff = clusterCutoff;

    return pNewInfoTable;
}

void SR_ReadPairInfoTableFree(SR_ReadPairInfoTable* pInfoTable)
{
    if (pInfoTable != NULL)
    {
        free(pInfoTable->chrIndex);
        free(pInfoTable->detectBound);
        free(pInfoTable->clusterBound);
        free(pInfoTable);
    }
}

void SR_ReadPairInfoTableSetCutOff(SR_ReadPairInfoTable* pInfoTable, SR_FragLenDstrb* pDstrb)
{
    SR_Bool foundDownDetectCutoff = FALSE;
    SR_Bool foundUpDetectCutoff = FALSE;
    SR_Bool foundDownClusterCutoff = FALSE;
    SR_Bool foundUpClusterCutoff = FALSE;

    for (unsigned int i = 0; i != pDstrb->size; ++i)
    {
        for (unsigned int j = 0; j != NUM_ALLOWED_HIST; ++j)
        {
            unsigned index = i * NUM_ALLOWED_HIST * 2;
            for (unsigned int k = 0; k != pDstrb->pHists[i].size[j]; ++k)
            {
                if (pDstrb->pHists[i].cdf[j][k] >= pInfoTable->detectCutoff && !foundDownDetectCutoff)
                {
                    pInfoTable->detectBound[index] = k;
                    foundDownDetectCutoff = TRUE;
                }

                if (pDstrb->pHists[i].cdf[j][k] >= pInfoTable->clusterCutoff && !foundDownClusterCutoff)
                {
                    pInfoTable->clusterBound[index] = pDstrb->pHists[i].fragLen[j][k];
                    foundDownClusterCutoff = TRUE;
                }

                if (foundDownClusterCutoff && foundDownDetectCutoff)
                    break;
            }

            for (int k = pDstrb->pHists[i].size[j] - 1; k != -1; --k)
            {
                if (pDstrb->pHists[i].cdf[j][k] <= 1 -  pInfoTable->detectCutoff && !foundUpDetectCutoff)
                {
                    pInfoTable->detectBound[index + 1] = k;
                    foundUpDetectCutoff = TRUE;
                }

                if (pDstrb->pHists[i].cdf[j][k] <= 1 - pInfoTable->clusterCutoff && !foundUpClusterCutoff)
                {
                    pInfoTable->clusterBound[index + 1] = pDstrb->pHists[i].fragLen[j][k];
                    foundUpClusterCutoff = TRUE;
                }

                if (foundUpClusterCutoff && foundUpDetectCutoff)
                    break;
            }
        }
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

    if (pUpAlgn->alignment.core.tid == pDownAlgn->alignment.core.tid)
        info.fragLen = abs(pUpAlgn->alignment.core.isize);
    else
    {
        info.upRefID = -(pUpAlgn->alignment.core.tid);
        info.eventType = SV_INTER_CHR_TRNSLCTN;
        return SR_OK;
    }

    int8_t histIndex = pDstrb->validModeMap[pairStats.pairMode];
    if (histIndex < 0)
    {
        if (pDstrb->numPairMode == 2)
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
    }
    else
    {
        // any reads with fragment length at the edge of the fragment length distribution will be kept as SV candidates
        uint32_t lowerCutoffIndex = pDstrb->pHists[info.readGrpID].cutoff[histIndex][DSTRB_LOWER_CUTOFF];
        uint32_t upperCutoffIndex = pDstrb->pHists[info.readGrpID].cutoff[histIndex][DSTRB_UPPER_CUTOFF];

        uint32_t lowerCutoff = pDstrb->pHists[info.readGrpID].fragLen[histIndex][lowerCutoffIndex];
        uint32_t upperCutoff = pDstrb->pHists[info.readGrpID].fragLen[histIndex][upperCutoffIndex];

        if (pairStats.fragLen < lowerCutoff)
        {
            SR_ReadPairLoadProb(&info, pDstrb, 0, lowerCutoffIndex);
            info.eventType = SV_INSERTION;
        }
        else if (pairStats.fragLen > upperCutoff)
        {
            SR_ReadPairLoadProb(&info, pDstrb, upperCutoffIndex, pDstrb->pHists[info.readGrpID].size[histIndex] - 1);
            info.eventType = SV_DELETION;
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
