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

#define NUM_TOTAL_PAIR_MODE 8

#define NUM_ALLOWED_PAIR_MODE 4

#define NUM_ALLOWED_HIST 2

#define INVALID_PAIR_MODE_SET_INDEX 2

typedef struct SR_FragLenHist
{
    void* rawHist[NUM_ALLOWED_HIST];

    uint32_t* fragLen[NUM_ALLOWED_HIST];

    double* cdf[NUM_ALLOWED_HIST];

    double mean[NUM_ALLOWED_HIST];

    double median[NUM_ALLOWED_HIST];

    double stdev[NUM_ALLOWED_HIST];

    uint32_t size[NUM_ALLOWED_HIST];

    uint64_t modeCount[NUM_ALLOWED_HIST + 1];

}SR_FragLenHist;

typedef struct SR_FragLenDstrb
{
    char** pReadGrpNames;

    void* pReadGrpHash;

    SR_FragLenHist* pHists;

    int8_t validModeMap[NUM_TOTAL_PAIR_MODE];
    
    uint8_t numPairMode;

    uint32_t size;

    uint32_t capacity;

    unsigned short minMQ;

    SR_Bool hasRG;

}SR_FragLenDstrb;

SR_FragLenDstrb* SR_FragLenDstrbAlloc(unsigned short minMQ, uint32_t capacity);

void SR_FragLenDstrbFree(SR_FragLenDstrb* pDstrb);

SR_Status SR_FragLenDstrbSetPairMode(SR_FragLenDstrb* pDstrb, const char* cmdArg);

SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader);

SR_Status SR_FragLenDstrbUpdate(SR_FragLenDstrb* pDstrb, const SR_BamPairStats* pPairStats);

void SR_FragLenDstrbFinalize(SR_FragLenDstrb* pDstrb);

void SR_FragLenDstrbWrite(const SR_FragLenDstrb* pDstrb, FILE* dstrbOutput);

#endif  /*SR_FRAGLENDSTRB_H*/
