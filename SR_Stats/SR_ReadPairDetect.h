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
#include "SR_ReadPairBuild.h"

typedef struct SR_ReadPairDetectPars
{
    char* workingDir;

    int workingRefID[2];

    int minNumClustered;

    int minEventLength;

}SR_ReadPairDetectPars;

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


void SR_ReadPairDetect(const SR_ReadPairDetectPars* pDetectPars);

void SR_LocalPairArrayRead(SR_LocalPairArray* pLongPairArray, FILE* input);

void SR_CrossPairArrayRead(SR_CrossPairArray* pCrossPairArray, FILE* input);

void SR_SpecialPairArrayRead(SR_SpecialPairArray* pSpecialPairArray, FILE* input);

void SR_SpecialPairTableReadID(SR_SpecialPairTable* pSpeicalPairTable, FILE* libInput);

// void SV_AssistArrayResize(SV_AssistArray* pAssistArray, unsigned int newSize);

/*  
void SR_ReadPairFindDel(SR_DelArray* pDelArray, SV_AssistArray* pAssistArray, const SR_LocalPairArray* pLongPairArray,
                        SR_Cluster* pDelCluster, const SR_LibInfoTable* pLibTable, const SR_ReadPairDetectPars* pPars);

void SR_DelEventGenotype(SR_DelArray* pDelArray, const SR_LibInfoTable* pLibTable);
*/


#endif  /*SR_READPAIRDETECT_H*/
