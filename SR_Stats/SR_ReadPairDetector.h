/*
 * =====================================================================================
 *
 *       Filename:  SR_ReadPairDetector.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/08/2011 11:44:06 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_READPAIRDETECTOR_H
#define  SR_READPAIRDETECTOR_H

#include "SR_FragLenDstrb.h"

#define SR_NUM_SV_TYPES 6

typedef enum
{
    SV_UNKNOWN            = 0,

    SV_DELETION           = 1,

    SV_INSERTION          = 2,

    SV_INVERSION          = 3,

    SV_TANDEM_DUP         = 4,

    SV_INTER_CHR_TRNSLCTN = 5,

    SV_MEI_INSERTION      = 6

}SV_EventType;

typedef struct SR_ReadPairInfo
{
    int32_t readGrpID;
    int32_t downRefID;

    union
    {
        int32_t upRefID;
        int32_t fragLen;
    };

    int32_t upPos;

    int32_t downPos;

    int32_t probPower:8, probNum:8, eventType:8, pairMode:8;

}SR_ReadPairInfo;

typedef struct SR_ReadPairInfoArray
{
    SR_ReadPairInfo* data;

    uint32_t size;

    uint32_t capacity;

}SR_ReadPairInfoArray;

typedef struct SR_ReadPairInfoTable
{
    SR_ReadPairInfoArray arrays[SR_NUM_SV_TYPES];

    uint32_t* chrIndex;
    
    uint32_t numChr;

    uint32_t numRG;

    uint32_t* detectBound;             // cutoff of the probability

    uint32_t* clusterBound;

    double* medianFragLen;

    double detectCutoff;

    double clusterCutoff;

}SR_ReadPairInfoTable;

SR_ReadPairInfoTable* SR_ReadPairInfoTableAlloc(uint32_t numChr, uint32_t numRG, double detectCutoff, double clusterCutoff);

void SR_ReadPairInfoTableFree(SR_ReadPairInfoTable* pInfoTable);

//=======================================================================
// function:
//      set the probability cutoff for all the histograms
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. cutoff: the fragment length probability cutoff
//========================================================================
void SR_ReadPairInfoTableSetCutOff(SR_ReadPairInfoTable* pInfoTable, SR_FragLenDstrb* pDstrb);

SR_Status SR_ReadPairInfoTableUpdate(SR_ReadPairInfoTable* pInfoTable, const SR_BamNode* pUpAlgn, const SR_BamNode* pDownAlgn, const SR_FragLenDstrb* pDstrb);

#endif  /*SR_READPAIRDETECTOR_H*/