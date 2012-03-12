/*
 * =====================================================================================
 *
 *       Filename:  SR_Cluster.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/25/2011 14:25:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_CLUSTER_H
#define  SR_CLUSTER_H

#include "SR_ReadPairAttrbt.h"

typedef struct SR_ClusterElmnt
{
    double mean[2];

    double std[2];

    double median[2];

    double min[2];

    double max[2];

    unsigned int numReadPair;

    unsigned int startIndex;

}SR_ClusterElmnt;

typedef struct SR_ClusterElmntArray
{
    SR_ClusterElmnt* data;

    unsigned int size;

    unsigned int capacity;

}SR_ClusterElmntArray;

typedef struct SR_Cluster
{
    const SR_ReadPairAttrbtArray* pAttrbtArray;

    SR_ClusterElmntArray* pElmntArray;

    int* pNext;

    int* pMap;
    
    int* pCount;

    unsigned int size;

    unsigned int capacity;

    int minReadPairNum;

    double minStd[2];

}SR_Cluster;

SR_Cluster* SR_ClusterAlloc(const SR_ReadPairAttrbtArray* pAttrbtArray, int minReadPairNum, double minStd[2]);

void SR_ClusterFree(SR_Cluster* pCluster);

void SR_ClusterReinit(SR_Cluster* pCluster, const SR_ReadPairAttrbtArray* pAttrbtArray);

void SR_ClusterMake(SR_Cluster* pCluster);

void SR_ClusterBuild(SR_Cluster* pCluster);

void SR_ClusterFinalize(SR_Cluster* pCluster);

int SR_ClusterMerge(SR_Cluster* pCluster, int dstIndex, int srcIndex);

#endif  /*SR_CLUSTER_H*/
