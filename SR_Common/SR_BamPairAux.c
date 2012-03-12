/*
 * =====================================================================================
 *
 *       Filename:  SR_BamPairAux.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/17/2011 09:43:55 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "SR_Error.h"
#include "SR_BamPairAux.h"
#include "SR_BamInStream.h"

static SR_Status SR_LoadMate(bam1_t* pMate, const bam1_t* pAlignment, bamFile pBamInput, bam_index_t* pBamIndex)
{
    // jump and read the first alignment in the given chromosome
    bam_iter_t pBamIter = bam_iter_query(pBamIndex, pAlignment->core.mtid, pAlignment->core.mpos, pAlignment->core.mpos + 1);

    int ret;
    SR_Status retStatus = SR_ERR;
    while ((ret = bam_iter_read(pBamInput, pBamIter, pMate)) >= 0)
    {
        if (strcmp(bam1_qname(pMate), bam1_qname(pAlignment)) == 0)
        {
            retStatus = SR_OK;
            break;
        }
    }

    bam_iter_destroy(pBamIter);
    return retStatus;
}


SR_FilterDataRP* SR_FilterDataRPAlloc(const SR_AnchorInfo* pAnchorInfo, uint32_t binLen)
{
    SR_FilterDataRP* pFilterData = (SR_FilterDataRP*) malloc(sizeof(SR_FilterDataRP));
    if (pFilterData == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a filter data object.\n");

    pFilterData->binLen = binLen;
    pFilterData->isFilled = FALSE;
    pFilterData->loadCross = FALSE;
    pFilterData->pAnchorInfo = pAnchorInfo;

    pFilterData->pBamInput = NULL;
    pFilterData->pBamIndex = NULL;

    pFilterData->pUpAlgn = bam_init1();
    pFilterData->pDownAlgn = bam_init1();

    return pFilterData;
}

void SR_FilterDataRPFree(SR_FilterDataRP* pFilterData)
{
    if (pFilterData != NULL)
    {
        if (pFilterData->pBamInput != NULL)
            bam_close(pFilterData->pBamInput);

        if (pFilterData->pBamIndex != NULL)
            bam_index_destroy(pFilterData->pBamIndex);

        bam_destroy1(pFilterData->pUpAlgn);
        bam_destroy1(pFilterData->pDownAlgn);

        free(pFilterData);
    }
}

void SR_FilterDataRPInit(SR_FilterDataRP* pFilterData, const char* bamInputFile)
{
    if(pFilterData->pBamInput != NULL)
        bam_close(pFilterData->pBamInput);

    if (pFilterData->pBamIndex != NULL)
        bam_index_destroy(pFilterData->pBamIndex);

    pFilterData->pBamInput = bam_open(bamInputFile, "r");
    if (pFilterData->pBamInput == NULL)
        SR_ErrQuit("ERROR: Cannot open the bam file: %s.\n", bamInputFile);

    pFilterData->pBamIndex = bam_index_load(bamInputFile);
    if (pFilterData->pBamIndex == NULL)
        SR_ErrQuit("ERROR: Cannot open the index file for this bam: %s.\n", bamInputFile);
}

SR_StreamCode SR_ReadPairFilter(const bam1_t* pAlignment, void* pFilterData, int32_t currRefID, int32_t currBinPos)
{
    // this is the fragment length distribution
    SR_FilterDataRP* pDataRP = (SR_FilterDataRP*) pFilterData;
    pDataRP->isFilled = FALSE;

    if ((pAlignment->core.flag & BAM_FPAIRED) == 0
        || strcmp(bam1_qname(pAlignment), "*") == 0
        || (pAlignment->core.flag & SR_READ_PAIR_FMASK) != 0)
    {
        return STREAM_PASS;
    }

    // if any of the mate of a read pair lands on an unused reference we filter them out
    if (pDataRP->pAnchorInfo->pLength[pAlignment->core.tid] < 0
        || pDataRP->pAnchorInfo->pLength[pAlignment->core.mtid] < 0)
    {
        return STREAM_PASS;
    }

    if (pAlignment->core.tid != pAlignment->core.mtid)
    {
        if (pDataRP->loadCross && pAlignment->core.tid < pAlignment->core.mtid)
        {
            SR_Status status = SR_LoadMate(pDataRP->pDownAlgn, pAlignment, pDataRP->pBamInput, pDataRP->pBamIndex);
            if (status == SR_OK)
            {
                bam_copy1(pDataRP->pUpAlgn, pAlignment);
                pDataRP->isFilled = TRUE;
                return STREAM_RETRUN;
            }
        }

        return STREAM_PASS;
    }

    if (currRefID != pAlignment->core.tid || currBinPos < 0)
    {
        currRefID = pAlignment->core.tid;
        currBinPos = pAlignment->core.pos;
    }

    SR_Bool shouldLoad = FALSE;
    if (pAlignment->core.pos <= pAlignment->core.mpos)
    {
        if (pAlignment->core.pos < currBinPos + pDataRP->binLen)
        {
            if (pAlignment->core.mpos >= currBinPos + 2 * pDataRP->binLen)
                shouldLoad = TRUE;
        }
        else if (pAlignment->core.pos >= currBinPos + 2 * pDataRP->binLen)
        {
            if (pAlignment->core.mpos >= pAlignment->core.pos + 2 * pDataRP->binLen)
                shouldLoad = TRUE;
        }
        else
        {
            if (pAlignment->core.mpos >= pAlignment->core.pos + 3 * pDataRP->binLen)
                shouldLoad = TRUE;
        }
    }
    else
    {
        if (pAlignment->core.pos < currBinPos + pDataRP->binLen)
        {
            if (pAlignment->core.mpos + pDataRP->binLen < currBinPos)
                return STREAM_PASS;
        }
        else if (pAlignment->core.pos >= 2 * currBinPos + pDataRP->binLen)
        {
            return STREAM_PASS;
        }
        else
        {
            if (pAlignment->core.mpos < currBinPos)
                return STREAM_PASS;
        }
    }

    if (shouldLoad)
    {
        SR_Status status = SR_LoadMate(pDataRP->pDownAlgn, pAlignment, pDataRP->pBamInput, pDataRP->pBamIndex);
        if (status == SR_OK)
        {
            bam_copy1(pDataRP->pUpAlgn, pAlignment);
            pDataRP->isFilled = TRUE;
            return STREAM_RETRUN;
        }

        return STREAM_PASS;
    }

    return STREAM_KEEP;
}

    /*
    // get the statistics of the read pair
    SR_PairStats pairStats;
    SR_Status status = SR_LoadPairStats(&pairStats, pAlignment);
    if (status == SR_ERR)
        return TRUE;

    // any reads do not have valid read group name will be filtered out
    int32_t readGrpIndex = 0;
    status = SR_FragLenDstrbGetRGIndex(&readGrpIndex, pDstrb, pairStats.RG);
    if (status == SR_ERR)
        return TRUE;
    
    // any reads aligned to different chromosome will be kept as SV candidates
    if (pAlignment->core.tid != pAlignment->core.mtid)
    {
        return FALSE;
    }

    // any reads aligned with improper pair mode (orientation) will be kept as SV candidates
    if (!SR_IsValidPairMode(pDstrb, pairStats.pairMode))
    {
        return FALSE;
    }

    // any reads with fragment length at the edge of the fragment length distribution will be kept as SV candidates
    uint32_t lowerCutoffIndex = pDstrb->pHists[readGrpIndex].cutoff[DSTRB_LOWER_CUTOFF];
    uint32_t upperCutoffIndex = pDstrb->pHists[readGrpIndex].cutoff[DSTRB_UPPER_CUTOFF];

    uint32_t lowerCutoff = pDstrb->pHists[readGrpIndex].fragLen[lowerCutoffIndex];
    uint32_t upperCutoff = pDstrb->pHists[readGrpIndex].fragLen[upperCutoffIndex];

    if (pairStats.fragLen < lowerCutoff || pairStats.fragLen > upperCutoff)
        return FALSE;

    // at last, those reads with valid pair mode and proper fragment length will be filtered out
    return TRUE;

    */


/*  

//====================================================================
// function:
//      check if a pair of read is normal
//
// args:
//      1. ppUpAlgn: a pointer of the pointer to a bam node object
//                   for the alignment with smaller coordinate
//      1. ppDownAlgn: a pointer of the pointer to a bam node object
//                   for the alignment with greater coordinate
//
// return:
//      if the pair of read is normal, return TRUE; else, return
//      FALSE
//=====================================================================
static inline SR_Bool SR_IsNormalPair(SR_BamNode** ppUpAlgn, SR_BamNode** ppDownAlgn, unsigned short minMQ)
{
    if (((*ppUpAlgn)->alignment.core.flag & BAM_FUNMAP) != 0
        || ((*ppDownAlgn)->alignment.core.flag & BAM_FUNMAP) != 0
        || (*ppUpAlgn)->alignment.core.qual < minMQ 
        || (*ppDownAlgn)->alignment.core.qual < minMQ)
    {
        return FALSE;
    }

    return TRUE;
}

*/
