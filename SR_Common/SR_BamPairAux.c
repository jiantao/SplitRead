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

#include "bam.h"
#include "SR_BamPairAux.h"


// default soft clipping tolerance
#define DEFAULT_SC_TOLERANCE 0.2

#define DEFAULT_MAX_MISMATCH_RATE 0.1

#define BAM_CMISMATCH 8

// alignment status
enum AlignmentStatus
{
    NONE_GOOD   = -1,       // neither a good anchor nor a good orphan candidate

    GOOD_ANCHOR    = 0,     // a good anchor candidate

    GOOD_ORPHAN    = 1,     // a good orphan candidate

    GOOD_SOFT      = 2,     // a good soft clipping candidate

    GOOD_MULTIPLE  = 3,     // a good multiple aligned candidate
};

static double SR_GetMismatchRate(const bam1_t* pAlignment)
{
    uint32_t* cigar = bam1_cigar(pAlignment);
    unsigned int numMismatch = 0;
    for (unsigned i = 0; i != pAlignment->core.n_cigar; ++i)
    {
        int type = (cigar[i] & BAM_CIGAR_MASK);
        if (type == BAM_CINS || type == BAM_CDEL || type == BAM_CSOFT_CLIP || type == BAM_CMISMATCH)
        {
            numMismatch += (cigar[i] >> BAM_CIGAR_SHIFT);
        }
    }

    return ((double) numMismatch / pAlignment->core.l_qseq);
}

static int SR_CheckAlignment(const bam1_t* pAlignment, double scTolerance, double maxMismatchRate, unsigned char minMQ)
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
        if (pAlignment->core.qual >= minMQ)
            return GOOD_ANCHOR;
        else if (maxMismatchRate > 0.0)
        {
            double mismatchRate = SR_GetMismatchRate(pAlignment);
            if (mismatchRate <= maxMismatchRate)
                return GOOD_MULTIPLE;
        }

        return NONE_GOOD;
    }
    else if (isHeadSC && isTailSC)
        return NONE_GOOD;
    else
        return GOOD_SOFT;
}

SR_AlgnType SR_GetAlignmentType(SR_BamNode** ppAlgnOne, SR_BamNode** ppAlgnTwo, double scTolerance, double maxMismatchRate, unsigned char minMQ)
{
    int firstType = SR_CheckAlignment(&((*ppAlgnOne)->alignment), scTolerance, maxMismatchRate, minMQ);
    int secondType = SR_CheckAlignment(&((*ppAlgnTwo)->alignment), scTolerance, maxMismatchRate, minMQ);

    if (firstType != GOOD_ANCHOR && secondType == GOOD_ANCHOR)
    {
        SR_SWAP(*ppAlgnOne, *ppAlgnTwo, SR_BamNode*);
        SR_SWAP(firstType, secondType, int);
    }

    if (firstType == GOOD_ANCHOR)
    {
        switch(secondType)
        {
            case GOOD_ORPHAN:
                return SR_UNIQUE_ORPHAN;
                break;
            case GOOD_SOFT:
                return SR_UNIQUE_SOFT;
                break;
            case GOOD_MULTIPLE:
                return SR_UNIQUE_MULTIPLE;
                break;
            case GOOD_ANCHOR:
                return SR_UNIQUE_NORMAL;
                break;
            default:
                return SR_OTHER_ALGN_TYPE;
                break;
        }
    }

    return SR_OTHER_ALGN_TYPE;
}

/* 
SR_Bool SR_IsUniqueOrphanPair(SR_BamNode** ppAnchor, SR_BamNode** ppOrphan, double scTolerance, unsigned char minMQ)
{
    int anchorStatus = SR_CheckAlignment(&((*ppAnchor)->alignment), scTolerance, minMQ);
    int orphanStatus = SR_CheckAlignment(&((*ppOrphan)->alignment), scTolerance, minMQ);

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

SR_Status SR_LoadUniquOrphanPairs(SR_BamInStream* pBamInStream, unsigned int threadID, double scTolerance, unsigned char minMQ)
{
    SR_BamNode* pAnchor = NULL;
    SR_BamNode* pOrphan = NULL;

    SR_Status readerStatus = SR_OK;
    SR_Status bufferStatus = SR_OK;
    while ((readerStatus = SR_BamInStreamLoadPair(&pAnchor, &pOrphan, pBamInStream)) == SR_OK)
    {
        if (SR_IsUniqueOrphanPair(&pAnchor, &pOrphan, scTolerance, minMQ))
        {
            SR_BamInStreamPush(pBamInStream, pAnchor, threadID);
            bufferStatus = SR_BamInStreamPush(pBamInStream, pOrphan, threadID);
        }
        else
        {
            SR_BamInStreamRecycle(pBamInStream, pAnchor);
            SR_BamInStreamRecycle(pBamInStream, pOrphan);
        }

        if (bufferStatus == SR_FULL)
            break;
    }

    return readerStatus;
}
*/

SR_Status SR_LoadAlgnPairs(SR_BamInStream* pBamInStream, SR_FragLenDstrb* pDstrb, unsigned int threadID, double scTolerance, double maxMismatchRate, unsigned char minMQ)
{
    SR_BamNode* pAlgnOne = NULL;
    SR_BamNode* pAlgnTwo = NULL;

    SR_Status readerStatus = SR_OK;
    SR_Status bufferStatus = SR_OK;
    while ((readerStatus = SR_BamInStreamLoadPair(&pAlgnOne, &pAlgnTwo, pBamInStream)) == SR_OK)
    {
        SR_AlgnType algnType = SR_GetAlignmentType(&pAlgnOne, &pAlgnTwo, scTolerance, maxMismatchRate, minMQ);
        if ((algnType == SR_UNIQUE_ORPHAN || algnType == SR_UNIQUE_SOFT || algnType == SR_UNIQUE_MULTIPLE) && pBamInStream->numThreads > 0)
        {
            SR_BamInStreamPush(pBamInStream, pAlgnOne, threadID);
            bufferStatus = SR_BamInStreamPush(pBamInStream, pAlgnTwo, threadID);

            SR_BamInStreamSetAlgnType(pBamInStream, threadID, algnType);
        }
        else
        {
            if (algnType == SR_UNIQUE_NORMAL && pDstrb != NULL)
            {
                SR_BamPairStats pairStats;
                SR_Status modeStatus = SR_LoadPairStats(&pairStats, pAlgnOne);
                if (modeStatus == SR_OK)
                    SR_FragLenDstrbUpdate(pDstrb, &pairStats);
            }

            SR_BamInStreamRecycle(pBamInStream, pAlgnOne);
            SR_BamInStreamRecycle(pBamInStream, pAlgnTwo);
        }

        if (bufferStatus == SR_FULL)
            break;
    }

    return readerStatus;
}
