/*
 * =====================================================================================
 *
 *       Filename:  SR_ReadPairData.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/28/2012 10:40:51 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "SR_Error.h"
#include "SR_ReadPairAttrbt.h"

#define DEFAULT_RP_ATTRB_CAPACITY 50

static int CompareAttrbt(const void* a, const void* b)
{
    const SR_ReadPairAttrbt* first = a;
    const SR_ReadPairAttrbt* second = b;

    if (first->firstAttribute < second->firstAttribute)
    {
        return -1;
    }
    else if (first->firstAttribute > second->firstAttribute)
    {
        return 1;
    }
    else
    {
        if (first->secondAttribute < second->secondAttribute)
            return -1;
        else if (first->secondAttribute > second->secondAttribute)
            return 1;
        else
        {
            if (first->readGrpID <= second->readGrpID)
                return -1;
            else
                return 1;
        }
    }
}


SR_ReadPairAttrbtArray* SR_ReadPairAttrbtArrayAlloc(uint32_t numReadGrp)
{
    SR_ReadPairAttrbtArray* pAttrbtArray = (SR_ReadPairAttrbtArray*) malloc(sizeof(SR_ReadPairAttrbtArray));
    if (pAttrbtArray == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a read pair attribute array object.\n");

    pAttrbtArray->data = (SR_ReadPairAttrbt*) malloc(DEFAULT_RP_ATTRB_CAPACITY * sizeof(SR_ReadPairAttrbt));
    if (pAttrbtArray->data == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the read pair attributes in the read pair attribute array object.\n");

    pAttrbtArray->pBoundaries = (double (*)[2]) malloc(sizeof(double) * 2 * numReadGrp);
    if (pAttrbtArray->pBoundaries == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the boundaries in the read pair attribute array object.\n");

    pAttrbtArray->numReadGrp = numReadGrp;
    pAttrbtArray->size = 0;
    pAttrbtArray->capapcity = DEFAULT_RP_ATTRB_CAPACITY;

    return pAttrbtArray;
}

void SR_ReadPairAttrbtArrayFree(SR_ReadPairAttrbtArray* pAttrbtArray)
{
    if (pAttrbtArray != NULL)
    {
        free(pAttrbtArray->data);
        free(pAttrbtArray->pBoundaries);
        free(pAttrbtArray);
    }
}


void SR_ReadPairAttrbtArrayReInit(SR_ReadPairAttrbtArray* pAttrbtArray, uint64_t newCapacity)
{
    pAttrbtArray->size = 0;

    if (newCapacity > pAttrbtArray->capapcity)
    {
        free(pAttrbtArray->data);
        pAttrbtArray->data = (SR_ReadPairAttrbt*) malloc(newCapacity * sizeof(SR_ReadPairAttrbt));
        if (pAttrbtArray->data == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the read pair attributes in the read pair attribute array object.\n");

        pAttrbtArray->capapcity = newCapacity;
    }
}

void SR_ReadPairMakeLocal(SR_ReadPairAttrbtArray* pAttrbtArray, const SR_LocalPairArray* pLocalPairArray, const SR_LibInfoTable* pLibTable, SV_ReadPairType readPairType)
{
    SR_ReadPairAttrbtArrayReInit(pAttrbtArray, pLocalPairArray->size);
    pAttrbtArray->readPairType = readPairType;
    pAttrbtArray->size = pLocalPairArray->size;
    pAttrbtArray->numReadGrp = pLibTable->size;

    for (unsigned int i = 0; i != pLocalPairArray->size; ++i)
    {
        pAttrbtArray->data[i].origIndex = i;
        pAttrbtArray->data[i].readGrpID = pLocalPairArray->data[i].readGrpID;

        double median = pLibTable->pLibInfo[pAttrbtArray->data[i].readGrpID].fraglenMedian;

        pAttrbtArray->data[i].firstAttribute =  pLocalPairArray->data[i].upPos + (double) pLocalPairArray->data[i].fragLen / 2;
        pAttrbtArray->data[i].secondAttribute = pLocalPairArray->data[i].fragLen - median;
    }

    qsort(pAttrbtArray->data, pAttrbtArray->size, sizeof(pAttrbtArray->data[0]), CompareAttrbt);

    // neighbourhood scale factor (borrowed from Spanner)
    double boundScale[2] = {1.25, 0.5};
    for (unsigned int i = 0; i != pAttrbtArray->numReadGrp; ++i)
    {
        pAttrbtArray->pBoundaries[i][0] = (double) pLibTable->pLibInfo[i].fraglenMedian * boundScale[0];
        pAttrbtArray->pBoundaries[i][1] = (double) (pLibTable->pLibInfo[i].fragLenHigh - pLibTable->pLibInfo[i].fragLenLow) * boundScale[1];
    }

    qsort(pAttrbtArray->data, pAttrbtArray->size, sizeof(pAttrbtArray->data[0]), CompareAttrbt);
}

void SR_ReadPairMakeCross(SR_ReadPairAttrbtArray* pAttrbtArray, const SR_LibInfoTable* pLibTable, const SR_CrossPairArray* pCrossPairArray)
{
    SR_ReadPairAttrbtArrayReInit(pAttrbtArray, pCrossPairArray->size);
    pAttrbtArray->readPairType = PT_CROSS;
    pAttrbtArray->size = pCrossPairArray->size;
    pAttrbtArray->numReadGrp = pLibTable->size;

    for (unsigned int i = 0; i != pCrossPairArray->size; ++i)
    {
        pAttrbtArray->data[i].origIndex = i;
        pAttrbtArray->data[i].readGrpID = pCrossPairArray->data[i].readGrpID;

        pAttrbtArray->data[i].firstAttribute = pCrossPairArray->data[i].upRefID * 1e10 + pCrossPairArray->data[i].upPos;
        pAttrbtArray->data[i].secondAttribute = pCrossPairArray->data[i].downRefID * 1e10 + pCrossPairArray->data[i].downPos;
    }

    qsort(pAttrbtArray->data, pAttrbtArray->size, sizeof(pAttrbtArray->data[0]), CompareAttrbt);

    double boundScale = 1.25;
    for (unsigned int i = 0; i != pAttrbtArray->numReadGrp; ++i)
    {
        pAttrbtArray->pBoundaries[i][0] = (double) (pLibTable->pLibInfo[i].fragLenHigh - pLibTable->pLibInfo[i].fragLenLow) * boundScale;
        pAttrbtArray->pBoundaries[i][1] = (double) (pLibTable->pLibInfo[i].fragLenHigh - pLibTable->pLibInfo[i].fragLenLow) * boundScale;
    }
}

void SR_ReadPairMakeInverted(SR_ReadPairAttrbtArray* pAttrbtArrays[2], const SR_LibInfoTable* pLibTable, const SR_LocalPairArray* pInvertedPairArray)
{
    SR_ReadPairAttrbtArrayReInit(pAttrbtArrays[0], pInvertedPairArray->subSize[0]);
    SR_ReadPairAttrbtArrayReInit(pAttrbtArrays[1], pInvertedPairArray->subSize[1]);

    pAttrbtArrays[0]->size = pInvertedPairArray->subSize[0];
    pAttrbtArrays[1]->size = pInvertedPairArray->subSize[1];

    pAttrbtArrays[0]->readPairType = PT_INVERTED3;
    pAttrbtArrays[1]->readPairType = PT_INVERTED5;

    pAttrbtArrays[0]->numReadGrp = pLibTable->size;
    pAttrbtArrays[1]->numReadGrp = pLibTable->size;

    for (unsigned int i = 0; i != pInvertedPairArray->size; ++i)
    {
        const SR_LocalPair* pInvertedPair = pInvertedPairArray->data + i;
        int arrayIndex = pInvertedPair->readPairType - PT_INVERTED3;

        pAttrbtArrays[arrayIndex]->data[i].origIndex = i;
        pAttrbtArrays[arrayIndex]->data[i].readGrpID = pInvertedPair->readGrpID;

        double median = pLibTable->pLibInfo[pInvertedPair->readGrpID].fraglenMedian;

        pAttrbtArrays[arrayIndex]->data[i].firstAttribute =  pInvertedPair->upPos + (double) pInvertedPair->fragLen / 2;
        pAttrbtArrays[arrayIndex]->data[i].secondAttribute = pInvertedPair->fragLen - median;
    }

    qsort(pAttrbtArrays[0]->data, pAttrbtArrays[0]->size, sizeof(pAttrbtArrays[0]->data[0]), CompareAttrbt);
    qsort(pAttrbtArrays[1]->data, pAttrbtArrays[1]->size, sizeof(pAttrbtArrays[1]->data[0]), CompareAttrbt);

    // neighbourhood scale factor (borrowed from Spanner)
    double boundScale[2] = {1.25, 0.5};
    for (unsigned int i = 0; i != pAttrbtArrays[0]->numReadGrp; ++i)
    {
        pAttrbtArrays[0]->pBoundaries[i][0] = (double) pLibTable->pLibInfo[i].fraglenMedian * boundScale[0] * 2.0;
        pAttrbtArrays[0]->pBoundaries[i][1] = (double) (pLibTable->pLibInfo[i].fragLenHigh - pLibTable->pLibInfo[i].fragLenLow) * boundScale[1];
    }

    memcpy(pAttrbtArrays[1]->pBoundaries, pAttrbtArrays[0]->pBoundaries, sizeof(double) * 2 * pAttrbtArrays[0]->numReadGrp);
}

void SR_ReadPairMakeSpecial(SR_ReadPairAttrbtArray* pAttrbtArrays[2], const SR_SpecialPairArray* pSpecialPairArray, const SR_LibInfoTable* pLibTable)
{
    SR_ReadPairAttrbtArrayReInit(pAttrbtArrays[0], pSpecialPairArray->subSize[0]);
    SR_ReadPairAttrbtArrayReInit(pAttrbtArrays[1], pSpecialPairArray->subSize[1]);

    pAttrbtArrays[0]->size = pSpecialPairArray->subSize[0];
    pAttrbtArrays[1]->size = pSpecialPairArray->subSize[1];

    pAttrbtArrays[0]->readPairType = PT_SPECIAL3;
    pAttrbtArrays[1]->readPairType = PT_SPECIAL5;

    for (unsigned int i = 0; i != pSpecialPairArray->size; ++i)
    {
        const SR_SpecialPair* pSpecialPair = pSpecialPairArray->data + i;
        int arrayIndex = pSpecialPair->readPairType - PT_SPECIAL3;
        double halfMedian = pLibTable->pLibInfo[pSpecialPair->readGrpID].fraglenMedian / 2.0;
        double halfMedians[2] = {halfMedian, -halfMedian};
        uint32_t pos[2] = {pSpecialPair->pos[0], pSpecialPair->end[0]};

        pAttrbtArrays[arrayIndex]->data[i].origIndex = i;
        pAttrbtArrays[arrayIndex]->data[i].readGrpID = pSpecialPair->readGrpID;

        pAttrbtArrays[arrayIndex]->data[i].firstAttribute = pos[arrayIndex] + halfMedians[arrayIndex];
        pAttrbtArrays[arrayIndex]->data[i].secondAttribute = 0;
    }

    qsort(pAttrbtArrays[0]->data, pAttrbtArrays[0]->size, sizeof(pAttrbtArrays[0]->data[0]), CompareAttrbt);
    qsort(pAttrbtArrays[1]->data, pAttrbtArrays[1]->size, sizeof(pAttrbtArrays[1]->data[0]), CompareAttrbt);

    double boundScale = 1.25;
    for (unsigned int i = 0; i != pAttrbtArrays[0]->numReadGrp; ++i)
    {
        pAttrbtArrays[0]->pBoundaries[i][0] = (double) (pLibTable->pLibInfo[i].fragLenHigh - pLibTable->pLibInfo[i].fragLenLow) * boundScale;
        pAttrbtArrays[1]->pBoundaries[i][1] = 1e-3;
    }

    memcpy(pAttrbtArrays[1]->pBoundaries, pAttrbtArrays[0]->pBoundaries, sizeof(double) * 2 * pAttrbtArrays[0]->numReadGrp);
}

