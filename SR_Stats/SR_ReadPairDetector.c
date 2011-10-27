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

#include "SR_Error.h"
#include "SR_ReadPairDetector.h"

static const unsigned int SR_READ_PAIR_FMASK = (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FMUNMAP);

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

SR_Bool SR_ReadPairFilter(SR_BamNode* pBamNode, const void* filterData)
{
    if ((pBamNode->alignment.core.flag & BAM_FPAIRED) == 0
        || strcmp(bam1_qname(&(pBamNode->alignment)), "*") == 0
        || (pBamNode->alignment.core.flag & SR_READ_PAIR_FMASK) != 0
        || (pBamNode->alignment.core.isize == 0))
    {
        return TRUE;
    }

    // get the statistics of the read pair
    SR_BamPairStats pairStats;
    SR_Status status = SR_BamPairStatsLoad(&pairStats, pBamNode);
    if (status == SR_ERR)
        return TRUE;

    // this is the fragment length distribution
    const SR_FragLenDstrb* pDstrb = (const SR_FragLenDstrb*) filterData;

    // any reads do not have valid read group name will be filtered out
    int32_t readGrpIndex = 0;
    status = SR_FragLenDstrbGetRGIndex(&readGrpIndex, pDstrb, pairStats.RG);
    if (status == SR_ERR)
        return TRUE;
    
    // any reads aligned to different chromosome will be kept as SV candidates
    if (pBamNode->alignment.core.tid != pBamNode->alignment.core.mtid)
    {
        return FALSE;
    }

    // any reads aligned with improper pair mode (orientation) will be kept as SV candidates
    int histIndex = pDstrb->validModeMap[pairStats.pairMode];
    if (histIndex < 0)
    {
        return FALSE;
    }

    // any reads with fragment length at the edge of the fragment length distribution will be kept as SV candidates
    uint32_t lowerCutoffIndex = pDstrb->pHists[readGrpIndex].cutoff[histIndex][DSTRB_LOWER_CUTOFF];
    uint32_t upperCutoffIndex = pDstrb->pHists[readGrpIndex].cutoff[histIndex][DSTRB_UPPER_CUTOFF];

    uint32_t lowerCutoff = pDstrb->pHists[readGrpIndex].fragLen[histIndex][lowerCutoffIndex];
    uint32_t upperCutoff = pDstrb->pHists[readGrpIndex].fragLen[histIndex][upperCutoffIndex];

    if (pairStats.fragLen < lowerCutoff || pairStats.fragLen > upperCutoff)
        return FALSE;

    // at last, those reads with valid pair mode and proper fragment length will be filtered out
    return TRUE;
}

SR_Status SR_ReadPairInfoLoad(SR_ReadPairInfo* pInfo, const SR_BamNode* pUpAlgn, const SR_BamNode* pDownAlgn, const SR_FragLenDstrb* pDstrb)
{
    pInfo->upPos = pUpAlgn->alignment.core.pos;
    pInfo->downPos = pDownAlgn->alignment.core.pos;

    // get the statistics of the read pair
    SR_BamPairStats pairStats;
    SR_BamPairStatsLoad(&pairStats, pUpAlgn);

    pInfo->pairMode = pairStats.pairMode;
    SR_FragLenDstrbGetRGIndex(&(pInfo->readGrpID), pDstrb, pairStats.RG);

    pInfo->eventType = SV_UNKNOWN;
    pInfo->probNum = 0;
    pInfo->probPower = 0;

    pInfo->upRefID = pUpAlgn->alignment.core.tid;

    if (pUpAlgn->alignment.core.tid == pDownAlgn->alignment.core.tid)
        pInfo->fragLen = abs(pUpAlgn->alignment.core.isize);
    else
    {
        pInfo->dowRefID = -(pDownAlgn->alignment.core.tid);
        pInfo->eventType |= SV_INTER_CHR_TRNSLCTN;
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
                pInfo->eventType = SV_UNKNOWN;
            }
            else if (eventFirst != SV_UNKNOWN)
            {
                pInfo->eventType |= eventFirst;
            }
            else
            {
                pInfo->eventType |= eventSecond;
            }
        }
        else
        {
            pInfo->eventType = SV_UNKNOWN;
        }
    }
    else
    {
        // any reads with fragment length at the edge of the fragment length distribution will be kept as SV candidates
        uint32_t lowerCutoffIndex = pDstrb->pHists[pInfo->readGrpID].cutoff[histIndex][DSTRB_LOWER_CUTOFF];
        uint32_t upperCutoffIndex = pDstrb->pHists[pInfo->readGrpID].cutoff[histIndex][DSTRB_UPPER_CUTOFF];

        uint32_t lowerCutoff = pDstrb->pHists[pInfo->readGrpID].fragLen[histIndex][lowerCutoffIndex];
        uint32_t upperCutoff = pDstrb->pHists[pInfo->readGrpID].fragLen[histIndex][upperCutoffIndex];

        if (pairStats.fragLen < lowerCutoff)
        {
            SR_ReadPairLoadProb(pInfo, pDstrb, 0, lowerCutoffIndex);
            pInfo->eventType |= SV_INSERTION;
        }
        else if (pairStats.fragLen > upperCutoff)
        {
            SR_ReadPairLoadProb(pInfo, pDstrb, upperCutoffIndex, pDstrb->pHists[pInfo->readGrpID].size[histIndex] - 1);
            pInfo->eventType |= SV_DELETION;
        }
        else
        {
            return SR_ERR;
        }
    }

    return SR_OK;
}
