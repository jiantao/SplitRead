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

    SV_TANDEM_DUP         = 3,

    SV_INVERSION3         = 4,

    SV_INVERSION5         = 5,

    SV_INTER_CHR_TRNSLCTN = 6,

    SV_MEI_ALU            = 7

}SV_EventType;

typedef struct SR_ReadPairInfo
{
    int32_t readGrpID;

    int32_t downRefID;

    int32_t upRefID;

    int32_t downPos;

    int32_t upPos;

    int32_t fragLen;

    int32_t probPower:8, probNum:8, eventType:8, pairMode:8;

}SR_ReadPairInfo;

typedef struct SR_ReadPairInfoArray
{
    SR_ReadPairInfo* data;

    uint32_t size;

    uint32_t capacity;

}SR_ReadPairInfoArray;

typedef struct SR_ReadPairAttrbt
{
    double firstAttribute;

    double secondAttribute;

    uint32_t origIdex;

    int32_t readGrpID;

}SR_ReadPairAttrbt;

typedef struct SR_ReadPairAttrbtArray
{
    SV_EventType eventType;

    SR_ReadPairAttrbt* data;

    uint32_t dataSize;

    uint32_t dataCap;

    double* bound;

    uint32_t boundSize;

    uint32_t boundCap;

}SR_ReadPairAttrbtArray;

typedef struct SR_ReadPairInfoTable
{
    SR_ReadPairInfoArray arrays[SR_NUM_SV_TYPES];

    uint32_t* chrIndex;
    
    uint32_t numChr;

    uint32_t numRG;

}SR_ReadPairInfoTable;

SR_ReadPairInfoTable* SR_ReadPairInfoTableAlloc(uint32_t numChr, uint32_t numRG, double cutoff);

void SR_ReadPairInfoTableFree(SR_ReadPairInfoTable* pInfoTable);

SR_ReadPairAttrbtArray* SR_ReadPairAttrbtArrayAlloc(uint32_t dataCap, uint32_t boundCap);

void SR_ReadPairAttrbtArrayFree(SR_ReadPairAttrbtArray* pAttrbtArray);

//=======================================================================
// function:
//      set the probability cutoff for all the histograms
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. cutoff: the fragment length probability cutoff
//========================================================================
void SR_ReadPairInfoTableSetCutOff(SR_ReadPairInfoTable* pInfoTable, const SR_FragLenDstrb* pDstrb);

SR_Status SR_ReadPairInfoTableUpdate(SR_ReadPairInfoTable* pInfoTable, const SR_BamNode* pUpAlgn, const SR_BamNode* pDownAlgn, const SR_FragLenDstrb* pDstrb);

void SR_ReadPairAttrbtArrayResize(SR_ReadPairAttrbtArray* pAttrbtArray, uint32_t newDataCap, uint32_t newBoundCap);

void SR_ReadPairMake(SR_ReadPairAttrbtArray* pAttrbtArray, const SR_FragLenDstrb* pDstrb, const SR_ReadPairInfoArray* pInfoArray, SV_EventType eventType);

#endif  /*SR_READPAIRDETECTOR_H*/
