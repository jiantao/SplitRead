/*
 * =====================================================================================
 *
 *       Filename:  SR_ReadPairData.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/28/2012 10:37:08 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_READPAIRDATA_H
#define  SR_READPAIRDATA_H

#include <stdint.h>

#include "SR_LibInfo.h"
#include "SR_ReadPairBuild.h"

typedef struct SR_ReadPairAttrbt
{
    double firstAttribute;

    double secondAttribute;

    uint64_t origIndex;

    int32_t readGrpID;

}SR_ReadPairAttrbt;

typedef struct SR_ReadPairAttrbtArray
{
    SR_ReadPairAttrbt* data;

    double (*pBoundaries)[2];

    uint64_t size;

    uint64_t capacity;

    uint32_t numReadGrp;

    SV_ReadPairType readPairType;

}SR_ReadPairAttrbtArray;

SR_ReadPairAttrbtArray* SR_ReadPairAttrbtArrayAlloc(uint32_t numReadGrp);

void SR_ReadPairAttrbtArrayFree(SR_ReadPairAttrbtArray* pAttrbtArray);

#define SR_ReadPairAttrbtArrayGetFirstBound(pAttrbtArray, i) ((pAttrbtArray)->pBoundaries[(pAttrbtArray)->data[(i)].readGrpID][0])

#define SR_ReadPairAttrbtArrayGetSecondBound(pAttrbtArray, i) ((pAttrbtArray)->pBoundaries[(pAttrbtArray)->data[(i)].readGrpID][1])

void SR_ReadPairAttrbtArrayReInit(SR_ReadPairAttrbtArray* pAttrbtArray, uint64_t newCapacity);

void SR_ReadPairMakeLocal(SR_ReadPairAttrbtArray* pAttrbtArray, const SR_LocalPairArray* pLocalPairArray, const SR_LibInfoTable* pLibTable, SV_ReadPairType readpairType);

void SR_ReadPairMakeCross(SR_ReadPairAttrbtArray* pAttrbtArray, const SR_LibInfoTable* pLibTable, const SR_CrossPairArray* pCrossPairArray);

void SR_ReadPairMakeInverted(SR_ReadPairAttrbtArray* pAttrbtArrays[2], const SR_LibInfoTable* pLibTable, const SR_LocalPairArray* pInvertedPairArray);

void SR_ReadPairMakeSpecial(SR_ReadPairAttrbtArray* pAttrbtArrays[2], const SR_SpecialPairArray* pSpecialArray, const SR_LibInfoTable* pLibTable);


#endif  /*SR_READPAIRDATA_H*/
