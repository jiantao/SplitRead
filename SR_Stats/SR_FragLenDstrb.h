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

#define NUM_PAIR_MODE 8

#define NUM_TOP_PAIR_MODE 2

typedef struct SR_PairModeBin
{
    SR_PairMode pairMode;

    uint64_t freq;

}SR_PairModeBin;

typedef struct SR_FragLenHist
{
    void* rawHist[NUM_PAIR_MODE];

    uint32_t* fragLen;

    double* cdf;

    double mean;

    double median;

    double stdev;

    uint32_t size;

    SR_PairModeBin modeCount[NUM_PAIR_MODE];

}SR_FragLenHist;

typedef struct SR_FragLenDstrb
{
    char** pReadGrpNames;

    void* pReadGrpHash;

    SR_FragLenHist* pHists;

    uint32_t size;

    uint32_t capacity;

    unsigned short minMQ;

    SR_Bool hasRG;

}SR_FragLenDstrb;

typedef struct SR_FragLenTable
{
    char** pReadGrpNames;

    void* pReadGrpHash;

    uint8_t* pairMode;

    uint64_t* totalFreq;

    double* mean;

    double* median;

    double* stdev;

    uint32_t* histIndex;

    uint32_t* fragLen;

    double*   cdf;

    uint32_t size;

    SR_Bool hasRG;

}SR_FragLenTable;

SR_FragLenDstrb* SR_FragLenDstrbAlloc(unsigned short minMQ, uint32_t capacity);

void SR_FragLenDstrbFree(SR_FragLenDstrb* pDstrb);

SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader);

SR_Status SR_FragLenDstrbUpdate(SR_FragLenDstrb* pDstrb, const SR_BamPairStats* pPairStats);

void SR_FragLenDstrbFinalize(SR_FragLenDstrb* pDstrb);

void SR_FragLenDstrbWrite(const SR_FragLenDstrb* pDstrb, FILE* dstrbOutput);



SR_FragLenTable* SR_FragLenTableAlloc(void);

void SR_FragLenTableFree(SR_FragLenTable* pFragLenTable);

void SR_FragLenTableRead(pFragLenTable* pFragLenTable, FILE* fragLenInput);

#endif  /*SR_FRAGLENDSTRB_H*/
