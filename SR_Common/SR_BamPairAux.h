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


#include "SR_Types.h"
#include "SR_BamMemPool.h"
#include "SR_BamInStream.h"
#include "SR_Utilities.h"


//===============================
// Type and constant definition
//===============================

// a map used to map the pair mode into its corresponding number
// negative value means invalid mode
static const int SR_PairModeMap[16] = { 
                                   -1, -1, 0, 1,

                                   -1, -1, 2, 3,

                                   4, 5, -1, -1,

                                   6, 7, -1, -1
                               };

// the object used to hold the basic statistics of a pair of alignments
typedef struct SR_BamPairStats
{
    const char* RG;           // name of the read group

    unsigned int fragLen;     // fragment length of the pair

    uint16_t pairMode;        // orientation mode of the pair

}SR_BamPairStats;


//======================
// Interface functions
//======================

//==================================================================
// function:
//      check if a pair of alignment is qualified unique-orphan 
//      pair
//
// args:
//      ppAnchor: a pointer to a pointer of bam node with anchor
//                alignment
//      ppOrphan: a pointer to a pointer of bam node with orphan
//                alignment
//      scTolerance: soft clipping tolerance
//
// return:
//      if they are qualified unique orphan pair return TRUE
//      else return FALSE
//==================================================================
SR_Bool SR_IsUniqueOrphanPair(SR_BamNode** ppAnchor, SR_BamNode** ppOrphan, double scTolerance);

//==================================================================
// function:
//      load a certain number of unique orphan pairs into a buffer
//      associated with a thread
//
// args:
//      1. pBamInStream: a pointer to a bam in stream object
//      2. threadID: the ID of the thread
//      3. scTolerance: soft clipping tolerance
//
// return:
//      status of the bam in stream. if we reach the end of a
//      chromosome, return SR_OUT_OF_RANGE; if we reach the end of
//      the file, return SR_EOF; if an error happens, return
//      SR_ERR; else return SR_OK
//==================================================================
SR_Status SR_LoadUniquOrphanPairs(SR_BamInStream* pBamInStream, unsigned int threadID, double scTolerance);

//==================================================================
// function:
//      get the orientation mode of a single read
//
// args:
//      1. pAlgn: a pointer to a bam alignment object
//
// return:
//      the mode of the single read
//==================================================================
inline SR_SingleMode SR_BamGetSingleMode(bam1_t* pAlgn)
{
    if ((pAlgn->core.flag & BAM_FREAD1) != 0)
    {
        if ((pAlgn->core.flag & BAM_FREVERSE) == 0)
            return SR_1F;
        else
            return SR_1R;
    }
    else
    {
        if ((pAlgn->core.flag & BAM_FREVERSE) == 0)
            return SR_2F;
        else
            return SR_2R;
    }
}

//==================================================================
// function:
//      get the pair orientation mode from two single read mode
//
// args:
//      1. upMode: mode of the read with samller coordinates
//      2. downMode: mode of the read with greater coordinates
//
// return:
//      the mode of the pair
//==================================================================
#define SR_BamGetPairMode(upMode, downMode) (((upMode) << 2) | (downMode))

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
inline SR_Bool SR_IsNormalPair(SR_BamNode** ppUpAlgn, SR_BamNode** ppDownAlgn, unsigned short minMQ)
{
    if (((*ppUpAlgn)->alignment.core.flag & BAM_FUNMAP) != 0
        || ((*ppDownAlgn)->alignment.core.flag & BAM_FUNMAP) != 0
        || (*ppUpAlgn)->alignment.core.qual < minMQ 
        || (*ppDownAlgn)->alignment.core.qual < minMQ)
    {
        return FALSE;
    }

    if ((*ppUpAlgn)->alignment.core.pos > (*ppDownAlgn)->alignment.core.pos)
    {
        SR_SWAP(*ppUpAlgn, *ppDownAlgn, SR_BamNode*);
    }

    return TRUE;
}

//=======================================================================
// function:
//      get the basic statistics from a pair of alignments
//
// args:
//      1. pPairStats: a poiter to the pair statistics object
//      2. ppUpAlgn: a pointer to a bam node object for the alignment 
//                   with smaller coordinate
//      3. ppDownAlgn: a pointer to a bam node object for the alignment 
//                     with greater coordinate
//
// return:
//      if an error occurs, return SR_ERR; else return SR_OK
//========================================================================
SR_Status SR_BamPairStatsLoad(SR_BamPairStats* pPairStats, SR_BamNode* pUpAlgn, SR_BamNode* pDownAlgn);

    
#endif  /*SR_BAMPAIRAUX_H*/
