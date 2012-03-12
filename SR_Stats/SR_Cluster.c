/*
 * =====================================================================================
 *
 *       Filename:  SR_Cluster.c
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  11/25/2011 23:01:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (),
 *        Company:
 *
 * =====================================================================================
 */

#include <math.h>

#include "khash.h"
#include "SR_Error.h"
#include "SR_Cluster.h"
#include "SR_Utilities.h"

#define DEFAULT_CLUSTER_SIZE 20

KHASH_MAP_INIT_INT(clusterMap, int);

static void SR_ClusterConnect(SR_Cluster* pCluster, unsigned int centerIndex, unsigned int memberIndex)
{
    if (centerIndex == memberIndex)
        return;

    int centerID = pCluster->pMap[centerIndex];
    int centerNext = pCluster->pNext[centerIndex];

    pCluster->pNext[centerIndex] = memberIndex;

    unsigned int i = memberIndex;
    for (; pCluster->pNext[i] != memberIndex; i = pCluster->pNext[i])
    {
        pCluster->pMap[i] = centerID;
    }

    pCluster->pNext[i] = centerNext;
    pCluster->pMap[i] = centerID;
}

static void SR_ClusterUpdateElmnt(SR_Cluster* pCluster, int clusterID, int attrbtID)
{
    SR_ClusterElmnt* pElmnt = SR_ARRAY_GET_PT(pCluster->pElmntArray, clusterID);
    SR_ReadPairAttrbt* pAttrbt = SR_ARRAY_GET_PT(pCluster->pAttrbtArray, attrbtID);

    ++pElmnt->numReadPair;

    pElmnt->mean[0] += pAttrbt->firstAttribute;
    pElmnt->mean[1] += pAttrbt->secondAttribute;

    pElmnt->std[0] += pAttrbt->firstAttribute * pAttrbt->firstAttribute;
    pElmnt->std[1] += pAttrbt->secondAttribute * pAttrbt->secondAttribute;

    pElmnt->min[0] = (pAttrbt->firstAttribute < pElmnt->min[0] ? pAttrbt->firstAttribute : pElmnt->min[0]);
    pElmnt->min[1] = (pAttrbt->secondAttribute < pElmnt->min[1] ? pAttrbt->secondAttribute : pElmnt->min[1]);

    pElmnt->max[0] = (pAttrbt->firstAttribute > pElmnt->max[0] ? pAttrbt->firstAttribute : pElmnt->max[0]);
    pElmnt->max[1] = (pAttrbt->secondAttribute > pElmnt->max[1] ? pAttrbt->secondAttribute : pElmnt->max[1]);
}

static int SR_ClusterAddElmnt(SR_Cluster* pCluster, int attrbtID)
{
    int clusterID = pCluster->pElmntArray->size;

    if (SR_ARRAY_IS_FULL(pCluster->pElmntArray))
    {
        SR_ARRAY_RESIZE(pCluster->pElmntArray, pCluster->pElmntArray->capacity * 2, SR_ClusterElmnt);
    }

    ++(pCluster->pElmntArray->size);

    SR_ClusterElmnt* pElmnt = SR_ARRAY_GET_PT(pCluster->pElmntArray, clusterID);
    SR_ReadPairAttrbt* pAttrbt = SR_ARRAY_GET_PT(pCluster->pAttrbtArray, attrbtID);

    pElmnt->numReadPair = 1;
    pElmnt->startIndex = attrbtID;

    pElmnt->mean[0] = pAttrbt->firstAttribute;
    pElmnt->mean[1] = pAttrbt->secondAttribute;

    pElmnt->std[0] = pAttrbt->firstAttribute * pAttrbt->firstAttribute;
    pElmnt->std[1] = pAttrbt->secondAttribute * pAttrbt->secondAttribute;

    pElmnt->min[0] = pAttrbt->firstAttribute;
    pElmnt->min[1] = pAttrbt->secondAttribute;

    pElmnt->max[0] = pAttrbt->firstAttribute;
    pElmnt->max[1] = pAttrbt->secondAttribute;

    return clusterID;
}

SR_Cluster* SR_ClusterAlloc(const SR_ReadPairAttrbtArray* pAttrbtArray, int minReadPairNum, double minStd[2])
{
    SR_Cluster* pCluster = (SR_Cluster*) malloc(sizeof(SR_Cluster));
    if (pCluster == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a cluster object.\n");

    SR_ARRAY_ALLOC(pCluster->pElmntArray, DEFAULT_CLUSTER_SIZE, SR_ClusterElmntArray, SR_ClusterElmnt);

    pCluster->pAttrbtArray = pAttrbtArray;
    pCluster->size = pAttrbtArray->size;
    pCluster->capacity = pAttrbtArray->size;

    pCluster->minReadPairNum = minReadPairNum;
    pCluster->minStd[0] = minStd[0];
    pCluster->minStd[1] = minStd[1];

    pCluster->pNext = (int*) malloc(pCluster->capacity * sizeof(int));
    if (pCluster->pNext == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of linked list array in a cluster object.\n");

    pCluster->pMap = (int*) malloc(pCluster->capacity * sizeof(int));
    if (pCluster->pMap == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the cluster map in a cluster object.\n");

    pCluster->pCount = (int*) malloc(pCluster->capacity * sizeof(int));
    if (pCluster->pCount == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of the neighbour count in a cluster object.\n");

    for (unsigned int i = 0; i != pCluster->size; ++i)
    {
        pCluster->pNext[i] = i;
        pCluster->pMap[i] = i;
        pCluster->pCount[i] = 0;
    }

    return pCluster;
}

void SR_ClusterFree(SR_Cluster* pCluster)
{
    if (pCluster != NULL)
    {
        free(pCluster->pNext);
        free(pCluster->pMap);
        free(pCluster->pCount);

        SR_ARRAY_FREE(pCluster->pElmntArray, TRUE);

        free(pCluster);
    }
}

void SR_ClusterReinit(SR_Cluster* pCluster, const SR_ReadPairAttrbtArray* pAttrbtArray)
{
    SR_ARRAY_RESET(pCluster->pElmntArray);

    pCluster->pAttrbtArray = pAttrbtArray;
    pCluster->size = pAttrbtArray->size;

    if (pAttrbtArray->size > pCluster->capacity)
    {
        free(pCluster->pNext);
        free(pCluster->pMap);

        pCluster->capacity = pAttrbtArray->size * 2;

        pCluster->pNext = (int*) malloc(pCluster->capacity * sizeof(int));
        if (pCluster->pNext == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of linked list array in a cluster object.\n");

        pCluster->pMap = (int*) malloc(pCluster->capacity * sizeof(int));
        if (pCluster->pMap == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the cluster map in a cluster object.\n");

        pCluster->pCount = (int*) malloc(pCluster->capacity * sizeof(int));
        if (pCluster->pCount == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of the neighbour count in a cluster object.\n");
    }

    for (unsigned int i = 0; i != pCluster->size; ++i)
    {
        pCluster->pNext[i] = i;
        pCluster->pMap[i] = i;
        pCluster->pCount[i] = 0;
    }
}

void SR_ClusterMake(SR_Cluster* pCluster)
{
    unsigned int lastLowIndex = 0;
    for (unsigned int i = 1; i != pCluster->size; ++i)
    {
        unsigned int j = lastLowIndex;

        double firstBound = SR_ReadPairAttrbtArrayGetFirstBound(pCluster->pAttrbtArray, i);
        double secondBound = SR_ReadPairAttrbtArrayGetSecondBound(pCluster->pAttrbtArray, i);

        while ((j < pCluster->size) && (pCluster->pAttrbtArray->data[i].firstAttribute - pCluster->pAttrbtArray->data[j].firstAttribute) > firstBound)
            ++j;

        lastLowIndex = j;

        while (fabs(pCluster->pAttrbtArray->data[i].firstAttribute - pCluster->pAttrbtArray->data[j].firstAttribute) <= firstBound)
        {
            if (fabs(pCluster->pAttrbtArray->data[i].secondAttribute - pCluster->pAttrbtArray->data[j].secondAttribute) <= secondBound)
                ++(pCluster->pCount[i]);

            ++j;

            if (j == pCluster->size)
                break;
        }
    }
}

void SR_ClusterBuild(SR_Cluster* pCluster)
{
    unsigned int lastLowIndex = 0;
    for (unsigned int i = 1; i != pCluster->size; ++i)
    {
        unsigned int j = lastLowIndex;

        double firstBound = SR_ReadPairAttrbtArrayGetFirstBound(pCluster->pAttrbtArray, i);
        double secondBound = SR_ReadPairAttrbtArrayGetSecondBound(pCluster->pAttrbtArray, i);

        while ((j < pCluster->size) && (pCluster->pAttrbtArray->data[i].firstAttribute - pCluster->pAttrbtArray->data[j].firstAttribute) > firstBound)
            ++j;

        lastLowIndex = j;

        unsigned int maxCount = 0;
        unsigned int centerIndex = 0;
        while (fabs(pCluster->pAttrbtArray->data[i].firstAttribute - pCluster->pAttrbtArray->data[j].firstAttribute) <= firstBound)
        {
            if (fabs(pCluster->pAttrbtArray->data[i].secondAttribute - pCluster->pAttrbtArray->data[j].secondAttribute) <= secondBound)
            {
                if (pCluster->pCount[j] > maxCount)
                {
                    centerIndex = j;
                    maxCount = pCluster->pCount[j];
                }
            }

            ++j;

            if (j == pCluster->size)
                break;
        }

        SR_ClusterConnect(pCluster, centerIndex, i);
    }
}


void SR_ClusterFinalize(SR_Cluster* pCluster)
{
    khash_t(clusterMap)* clusterHash = kh_init(clusterMap);
    kh_resize(clusterMap, clusterHash, DEFAULT_CLUSTER_SIZE);

    int khRet = 0;
    khiter_t khIter = 0;
    int clusterIndex = 0;

    for (unsigned int i = 0; i != pCluster->size; ++i)
    {
        khIter = kh_put(clusterMap, clusterHash, pCluster->pMap[i], &khRet);
        if (khRet == 0)
        {
            clusterIndex = kh_value(clusterHash, khIter);
            SR_ClusterUpdateElmnt(pCluster, clusterIndex, i);
        }
        else
        {
            clusterIndex = SR_ClusterAddElmnt(pCluster, i);
            kh_value(clusterHash, khIter) = clusterIndex;
        }

        pCluster->pMap[i] = clusterIndex;
    }

    kh_destroy(clusterMap, clusterHash);

    for (unsigned int i = 0; i != pCluster->pElmntArray->size; ++i)
    {
        SR_ClusterElmnt* pElmnt = SR_ARRAY_GET_PT(pCluster->pElmntArray, i);

        for (unsigned int j = 0; j != 2; ++j)
        {
            pElmnt->mean[j] /= pElmnt->numReadPair;
            pElmnt->std[j] = sqrt(pElmnt->std[j] / pElmnt->numReadPair - (pElmnt->mean[j] * pElmnt->mean[j]));
        }
    }
}

int SR_ClusterClean(SR_Cluster* pCluster)
{
    int oldReadPairNum = pCluster->pElmntArray->size;

    // lazy deletion
    // the deleted elements' number of read pairs will be set to zero
    for (unsigned int i = 0; i != pCluster->pElmntArray->size; ++i)
    {
        if (pCluster->pElmntArray->data[i].numReadPair <= pCluster->minReadPairNum
            || pCluster->pElmntArray->data[i].std[0] <= pCluster->minStd[0]
            || pCluster->pElmntArray->data[i].std[1] <= pCluster->minStd[1])
        {
            pCluster->pElmntArray->data[i].numReadPair = 0;
            --(pCluster->pElmntArray->size);
        }
    }

    return oldReadPairNum;
}

int SR_ClusterMerge(SR_Cluster* pCluster, int dstIndex, int srcIndex)
{
    SR_ClusterElmnt* pSrcElmnt = pCluster->pElmntArray->data + srcIndex;
    SR_ClusterElmnt* pDstElmnt = pCluster->pElmntArray->data + dstIndex;

    SR_ClusterConnect(pCluster, pDstElmnt->startIndex, pSrcElmnt->startIndex);

    for (unsigned int i = 0; i != 2; ++i)
    {
        pDstElmnt->min[i] = (pDstElmnt->min[i] < pSrcElmnt->min[i] ? pDstElmnt->min[i] : pSrcElmnt->min[i]);
        pDstElmnt->max[i] = (pDstElmnt->max[i] > pSrcElmnt->max[i] ? pDstElmnt->max[i] : pSrcElmnt->max[i]);

        pDstElmnt->std[i] = (pDstElmnt->std[i] * pDstElmnt->std[i] + pDstElmnt->mean[i] + pDstElmnt->mean[i]) * pDstElmnt->numReadPair;
        pSrcElmnt->std[i] = (pSrcElmnt->std[i] * pSrcElmnt->std[i] + pSrcElmnt->mean[i] + pSrcElmnt->mean[i]) * pSrcElmnt->numReadPair;

        pDstElmnt->std[i] += pSrcElmnt->std[i];

        pDstElmnt->mean[i] = (pDstElmnt->mean[i] * pDstElmnt->numReadPair + pSrcElmnt->mean[i] * pSrcElmnt->numReadPair) / (pDstElmnt->numReadPair + pSrcElmnt->numReadPair);
    }

    pDstElmnt->numReadPair += pSrcElmnt->numReadPair;
    pSrcElmnt->numReadPair = 0;
    --(pCluster->size);

    pDstElmnt->std[0] = sqrt(pDstElmnt->std[0] / pDstElmnt->numReadPair - pDstElmnt->mean[0] * pDstElmnt->mean[0]);
    pDstElmnt->std[1] = sqrt(pDstElmnt->std[1] / pDstElmnt->numReadPair - pDstElmnt->mean[1] * pDstElmnt->mean[1]);

    return pDstElmnt->numReadPair;
}
