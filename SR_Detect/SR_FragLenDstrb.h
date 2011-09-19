/*
 * =====================================================================================
 *
 *       Filename:  SR_FragLenDstrb.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/07/2011 03:06:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef  SR_FRAGLENDSTRB_H
#define  SR_FRAGLENDSTRB_H

#include "bam.h"
#include "SR_Types.h"
#include "SR_BamHeader.h"
#include "SR_BamPairAux.h"

typedef struct SR_FragLenHist
{
    uint64_t* bin;

    uint32_t size;

    uint32_t lowerBound;

    uint32_t min;

    uint32_t max;

    uint32_t mode;

    uint64_t total;

}SR_FragLenHist;

typedef struct SR_FragLenDstrb
{
    char** pReadGrpNames;

    void* pReadGrpHash;

    SR_FragLenHist* pHists;

    uint32_t size;

    uint32_t capacity;

    uint64_t pairModeCount[8];

    unsigned short minMQ;

    SR_Bool hasRG;

}SR_FragLenDstrb;

uint32_t SR_FragLenHistGetMean(const SR_FragLenHist* pHist);

uint32_t SR_FragLenHistGetMedian(const SR_FragLenHist* pHist);

SR_FragLenDstrb* SR_FragLenDstrbAlloc(unsigned short minMQ, uint32_t capacity);

void SR_FragLenDstrbFree(SR_FragLenDstrb* pDstrb);

SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader);

SR_Status SR_FragLenDstrbUpdate(SR_FragLenDstrb* pDstrb, const SR_BamPairStats* pPairStats);

#endif  /*SR_FRAGLENDSTRB_H*/
