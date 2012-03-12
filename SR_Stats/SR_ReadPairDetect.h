/*
 * =====================================================================================
 *
 *       Filename:  SR_ReadPairDetect.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/23/2012 01:38:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_READPAIRDETECT_H
#define  SR_READPAIRDETECT_H

#include "SR_Cluster.h"

typedef struct SR_DetectControlPars
{
    int minNumClustered;

    int minEventLength;

}SR_DetectControlPars;

typedef struct SR_DelEvent
{
    int32_t clusterID;

    int32_t refID;
    uint32_t pos;
    uint32_t end;
    uint32_t length;

    int pos5[3];
    int pos3[3];

    int CIpos[2];
    int CIend[2];
    int CIlen[2];

    unsigned char quality;
    unsigned char mapQ5;
    unsigned char mapQ3;

}SR_DelEvent;

typedef struct SR_DelArray
{
    SR_DelEvent* data;

    unsigned int size;

    unsigned int capacity;

}SR_DelArray;

typedef struct SV_AssistArray
{
    int* pFragLenDiff;

    int* pMapQ5;

    int* pMapQ3;

    unsigned int size;

    unsigned int capacity;

}SV_AssistArray;

SV_AssistArray* SV_AssistArrayAlloc(void);

void SV_AssistArrayFree(SV_AssistArray* pAssistArray);

void SV_AssistArrayResize(SV_AssistArray* pAssistArray, unsigned int newSize);

void SR_ReadPairFindDel(SR_DelArray* pDelArray, SV_AssistArray* pAssistArray, const SR_LocalPairArray* pLongPairArray,
                        SR_Cluster* pDelCluster, const SR_LibInfoTable* pLibTable, const SR_DetectControlPars* pPars);


void SR_DelEventGenotype(SR_DelArray* pDelArray, const SR_LibInfoTable* pLibTable);


#endif  /*SR_READPAIRDETECT_H*/
