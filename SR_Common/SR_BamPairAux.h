/*
 * =====================================================================================
 *
 *       Filename:  SR_BamPairAux.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/17/2011 09:43:53 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_BAMPAIRAUX_H
#define  SR_BAMPAIRAUX_H


#include "bam.h"
#include "SR_Types.h"
#include "SR_Utilities.h"
#include "SR_LibInfo.h"


//===============================
// Type and constant definition
//===============================

static const unsigned int SR_UNIQUE_ORPHAN_FMASK = (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);

static const unsigned int SR_NORMAL_FMASK = (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FMUNMAP);

static const unsigned int SR_READ_PAIR_FMASK = (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FMUNMAP);

typedef struct SR_FilterDataRP
{
    bam1_t* pUpAlgn;

    bam1_t* pDownAlgn;

    bamFile pBamInput;

    bam_index_t* pBamIndex;

    const SR_AnchorInfo* pAnchorInfo;

    uint32_t binLen;

    SR_Bool loadCross;

    SR_Bool isFilled;

}SR_FilterDataRP;


SR_FilterDataRP* SR_FilterDataRPAlloc(const SR_AnchorInfo* pAnchorInfo, uint32_t binLen);

void SR_FilterDataRPFree(SR_FilterDataRP* pFilterData);

//======================
// Interface functions
//======================

static inline SR_StreamCode SR_CommonFilter(const bam1_t* pAlignment, void* pFilterData, int32_t currRefID, int32_t currBinPos)
{
    if ((pAlignment->core.flag & BAM_FPAIRED) == 0
        || strcmp(bam1_qname(pAlignment), "*") == 0
        || (pAlignment->core.flag & SR_UNIQUE_ORPHAN_FMASK) != 0
        || (pAlignment->core.flag | (BAM_FUNMAP | BAM_FMUNMAP)) == pAlignment->core.flag)
    {
        return STREAM_PASS;
    }
    
    return STREAM_KEEP;
}

static inline SR_StreamCode SR_NormalFilter(const bam1_t* pAlignment, void* pFilterData, int32_t currRefID, int32_t currBinPos)
{
    if ((pAlignment->core.flag & BAM_FPAIRED) == 0
        || strcmp(bam1_qname(pAlignment), "*") == 0
        || (pAlignment->core.flag & SR_NORMAL_FMASK) != 0
        || (pAlignment->core.isize == 0))
    {
        return STREAM_PASS;
    }

    // any reads aligned to different chromosome will be kept as SV candidates
    if (pAlignment->core.tid != pAlignment->core.mtid)
        return STREAM_PASS;

    return STREAM_KEEP;
}

void SR_FilterDataRPInit(SR_FilterDataRP* pFilterData, const char* bamInputFile);

#define SR_FilterDataRPTurnOffCross(pFilterData) (pFilterData)->loadCross = FALSE

#define SR_FilterDataRPTurnOnCross(pFilterData) (pFilterData)->loadCross = TRUE

SR_StreamCode SR_ReadPairFilter(const bam1_t* pAlignment, void* pFilterData, int32_t cuuRefID, int32_t currBinPos);


    
#endif  /*SR_BAMPAIRAUX_H*/
